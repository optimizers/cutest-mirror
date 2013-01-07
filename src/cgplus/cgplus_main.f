C     ( Last modified on 2 Jan 2013 at 15:10:00 )
      PROGRAM          CGPMA
C
C  CG+ test driver for problems derived from SIF files.
C
C  Nick Gould, for CGT Productions.
C  July 2004
C  Revised for CUTEst, January 2013
C
      INTEGER          out  , N , INPUT , MAXIT, status 
      INTEGER          LP    , MP, I , METHOD, ITER  , NFUN, IREST
      INTEGER          IFLAG , INSPEC, IPRINT( 2 )
      INTEGER :: io_buffer = 11
      DOUBLE PRECISION F, EPS, GNORM , BIGINF, ZERO, TLEV
      LOGICAL          BOUNDS, FINISH
      PARAMETER      ( out  = 6 )
      PARAMETER      ( INPUT = 55, INSPEC = 56 )
      PARAMETER      ( BIGINF = 9.0D+19, ZERO = 0.0D0 )
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : ) :: X, G, D, GOLD, W
      CHARACTER ( LEN = 10 ) :: PNAME, SPCDAT
      CHARACTER ( LEN = 10 ), ALLOCATABLE, DIMENSION( : )  :: XNAMES
      COMMON / CGDD /  MP, LP
      COMMON / RUNINF / ITER, NFUN
      DOUBLE PRECISION :: CPU( 2 ), CALLS( 4 )
C     
C  Open the Spec file for the method.
C
      SPCDAT = 'CGPLUS.SPC'
      OPEN ( INSPEC, FILE = SPCDAT, FORM = 'FORMATTED',
     *      STATUS = 'OLD' )
      REWIND INSPEC
C
C  Read input Spec data.
C
C     IPRINT(1): specifies the frequency of output
C     IPRINT(2): specifies the amount of output
C     METHOD   : method used (1=Fletcher-Reeves,2=Polak-Ribiere,3=P-R+)
C     IREST    : no restart (0) or restart every n iterations (1)
C     MAXIT    : maximum number of iterations
C     EPS      : the required norm of the gradient
C
      READ ( INSPEC, 1000 ) IPRINT( 1 ), IPRINT( 2 ), METHOD, IREST,
     *                      MAXIT, EPS
C
C  Close input file.
C
      CLOSE ( INSPEC )
C
C  Open the relevant file.
C
      OPEN ( INPUT, FILE = 'OUTSDIF.d', FORM = 'FORMATTED',
     *       STATUS = 'OLD' )
C
C  Check to see if there is sufficient room
C
      CALL CUTEST_udimen( status, INPUT, N )
      IF ( status /= 0 ) GO TO 910

      ALLOCATE( X( n ), G( n ), D( n ), GOLD( n ), W( n ), XNAMES( n ), 
     *          STAT = status )
      IF ( status /= 0 ) GO TO 990
C
C  Set up SIF data.
C
      CALL CUTEST_usetup( status, INPUT, out, io_buffer, N, X, W, GOLD )
      IF ( status /= 0 ) GO TO 910
C
C  Obtain variable names.
C
      CALL CUTEST_unames( status, N, PNAME, XNAMES )
      IF ( status /= 0 ) GO TO 910
C
C  Set up algorithmic input data.
C
      BOUNDS  = .FALSE.
      DO 10 I = 1, N
         IF ( W( I ) .GT. - BIGINF .OR. GOLD( I ) .LT. BIGINF )
     *      BOUNDS = .TRUE.
   10 CONTINUE
      IF ( BOUNDS ) WRITE( out, 2030 )
      LP     = out
      MP     = out
      ITER  = - 1
      IFLAG  = 0
      FINISH = .FALSE.
C
C  Optimization loop
C
   20 CONTINUE
         ITER = ITER + 1
C
C  Evaluate the function and gradient
C
         CALL CUTEST_uofg( status, N, X, F, G, .TRUE. )
         IF ( status /= 0 ) GO TO 910
C
C  Call the optimizer.
C
   30    CONTINUE
         CALL CGFAM( N, X, F, G, D, GOLD, IPRINT, EPS, W,
     *               IFLAG, IREST, METHOD, FINISH )
C
C  Test for termination
C
        IF ( IFLAG .LE. 0 .OR. ITER .GT. MAXIT ) GO TO 50
        IF ( IFLAG .EQ. 1 ) GO TO 20
        IF ( IFLAG .EQ. 2 ) THEN
C
C Termination Test.  The user may replace it by some other test. However, 
C the parameter 'FINISH' must be set to 'TRUE' when the test is satisfied.
C
           TLEV = EPS * ( 1.0D+0 + ABS( F ) )
           DO 40 I = 1, N
              IF( ABS( G( I ) ) .GT. TLEV ) GO TO 30
  40       CONTINUE
           FINISH = .TRUE.
           GO TO 30
        ENDIF

   50 CONTINUE
C
C  Terminal exit.
C
      CALL CUTEST_ureport( status, CALLS, CPU )
      IF ( status /= 0 ) GO TO 910
      GNORM    = ZERO
      DO 110 I  = 1, N
         GNORM = MAX( GNORM, ABS( G( I ) ) )
  110 CONTINUE
      WRITE ( out, 2010 ) F, GNORM
C      DO 120 I = 1, N
C         WRITE( out, 2020 ) XNAMES( I ), X( I ), G( I )
C  120 CONTINUE
      WRITE ( out, 2000 ) PNAME, N, INT( CALLS(1) ), INT( CALLS(2) ),
     *                     IFLAG, F, CPU(1), CPU(2) 
      CLOSE( INPUT  )
      CALL CUTEST_uterminate( status )
      STOP

  910 CONTINUE
      WRITE( out, "( ' CUTEst error, status = ', i0, ', stopping' )") 
     *   status
      STOP

  990 CONTINUE
      WRITE( out, "( ' Allocation error, status = ', I0 )" ) status
      STOP
C
C  Non-executable statements.
C
 1000 FORMAT( 5( I10, / ), D10.3 )
 2000 FORMAT( /, 24('*'), ' CUTEr statistics ', 24('*') //
     *    ,' Code used               :  CG+',     /
     *    ,' Problem                 :  ', A10,    /
     *    ,' # variables             =      ', I10 /
     *    ,' # objective functions   =      ', I10 /
     *    ,' # objective gradients   =      ', I10 / 
     *     ' Exit code               =      ', I10 /
     *    ,' Final f                 = ', E15.7 /
     *    ,' Set up time             =      ', 0P, F10.2, ' seconds' /
     *     ' Solve time              =      ', 0P, F10.2, ' seconds' //
     *     66('*') / )
 2010 FORMAT( ' Final objective function value  = ', 1P, D12.4, 
     *        /, ' Final norm of gradient          = ', 1P, D12.4,
     *        //, '                 X         G ' )
 2020 FORMAT(  1X, A10, 1P, 2D12.4 )
 2030 FORMAT(  /, ' ** Warning from CGPMA. The problem as stated',
     *            ' includes simple bounds. ', /,
     *            '    These bounds will be ignored. ' )
      END
