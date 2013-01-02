C     ( Last modified on 2 Jan 2013 at 15:10:00 )
      PROGRAM          CGDMA
C
C  CG_DESCENT test driver for problems derived from SIF files.
C
C  Nick Gould, for CGT Productions.
C  July 2004
C  Revised for CUTEst, January 2013
C
C  A number of updates to the structure of the SPEC file
C  to accomodate version 1.4 of CG_DESCENT.
C  Dominique Orban, July 2007
C
      IMPLICIT NONE
      INTEGER          IOUT  , N , M , INPUT , MSAVE, IOUTCP
      INTEGER          LP    , MP, I , METHOD, ITER  , NF, NG
      INTEGER          STAT  , ICALL , INSPEC, IPRINT( 2 )
      INTEGER :: io_buffer = 11
      DOUBLE PRECISION F, TOL, GNORM , BIGINF, ZERO,
     *                 STPMIN, STPMAX, TLEV
      LOGICAL          BOUNDS
      DOUBLE PRECISION DELTA, SIGMA, EPSILON, THETA, GAMMA,STOPFA,
     *                 RHO, ETA, PSI0, PSI1, PSI2, QUADCU, RSTRTF, 
     *                 MAXITF, FEPS, AWLFFCT, QDECAY
      INTEGER          NXPAND, NSECNT
      LOGICAL          QUADST, PRNTLV, PRNTFI, STRULE, ERULE, AWOLFE, 
     *                 STEP, PRTRUL, DEBUG
      PARAMETER      ( IOUT  = 6 )
      PARAMETER      ( INPUT = 55, INSPEC = 56, IOUTCP = 57 )
      PARAMETER      ( BIGINF = 9.0D+19, ZERO = 0.0D0 )
 
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : ) :: X, G, D, XTEMP, GTEMP
      CHARACTER ( LEN = 10 ) :: PNAME
      CHARACTER ( LEN = 10 ), ALLOCATABLE, DIMENSION( : )  :: XNAMES

 
      CHARACTER * 17   CGPARM
      PARAMETER      ( CGPARM = 'cg_descent_f.parm' )
      CHARACTER * 15   SPCDAT
      REAL             CPU( 2 ), CALLS( 4 )
      EXTERNAL         EVALF, EVALG
C     
C  Open the Spec file for the method.
C
      SPCDAT = 'CG_DESCENT.SPC'
      OPEN ( INSPEC, FILE = SPCDAT, FORM = 'FORMATTED',
     *      STATUS = 'OLD' )
      REWIND INSPEC
C
C  Read input Spec data.
C
      READ ( INSPEC, 1000 ) TOL, GNORM, DELTA, SIGMA, EPSILON, GAMMA,
     *                      RHO, ETA,
     *                      PSI0, PSI1, PSI2, QUADCU, STOPFA, AWLFFCT,
     *                      RSTRTF, MAXITF, FEPS, QDECAY, NXPAND,
     *                      NSECNT, PRTRUL, QUADST, PRNTLV, PRNTFI,
     *                      STRULE, AWOLFE, STEP, DEBUG
C
C  Close input file.
C
      CLOSE ( INSPEC )
C
C  Create the required data input file
C
      OPEN ( IOUTCP, FILE = CGPARM, FORM = 'FORMATTED',
     *       STATUS = 'UNKNOWN' )
      REWIND IOUTCP
      WRITE( IOUTCP, 1001 ) DELTA, SIGMA, EPSILON, GAMMA, RHO, ETA,
     *                      PSI0, PSI1, PSI2, QUADCU, STOPFA, AWLFFCT,
     *                      RSTRTF, MAXITF, FEPS, QDECAY, NXPAND,
     *                      NSECNT, PRTRUL, QUADST, PRNTLV, PRNTFI,
     *                      STRULE, AWOLFE, STEP, DEBUG
      CLOSE ( IOUTCP )
C
C  Open the relevant file.
C
      OPEN ( INPUT, FILE = 'OUTSDIF.d', FORM = 'FORMATTED',
     *       STATUS = 'OLD' )
C
C  Check to see if there is sufficient room
C
      CALL CUTEST_udimen( INPUT, status, N )
      IF ( status /= 0 ) GO TO 910

      ALLOCATE( X( n ), G( n ), D( n ), XTEMP( n ), GTEMP( n ), 
     *          XNAMES( n ), STAT = status )
      IF ( status /= 0 ) GO TO 990
C
C  Set up SIF data.
C
      CALL CUTEST_usetup( status, INPUT, IOUT, N, X, XTEMP, GTEMP )
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
         IF ( XTEMP( I ) .GT. - BIGINF .OR. GTEMP( I ) .LT. BIGINF )
     *      BOUNDS = .TRUE.
   10 CONTINUE
      IF ( BOUNDS ) WRITE( IOUT, 2030 )
      LP     = IOUT
      MP     = IOUT
C
C  Set up initial step length if requested
C
      IF( STEP ) THEN
         IF( GNORM .LE. 0.0D+0 ) THEN
             CALL UGR( N, X, G )
             GNORM = 0.0D+0
             DO 11 I = 1, N
                GNORM = MAX( GNORM, DABS(G(I)) )
   11        CONTINUE
         ENDIF
      ENDIF
C
C  Call the optimizer.
C
      CALL CG_DESCENT( TOL, X, N, EVALF, EVALG, STAT, GNORM, F,
     &                 ITER, NF, NG, D, G, XTEMP, GTEMP )
C
C  Terminal exit.
C
      CALL CUTEST_ureport( status, CALLS, CPU )
      IF ( status /= 0 ) GO TO 910
      WRITE ( IOUT, 2010 ) F, GNORM
C      DO 120 I = 1, N
C         WRITE( IOUT, 2020 ) XNAMES( I ), X( I ), G( I )
C  120 CONTINUE
      WRITE ( IOUT, 2000 ) PNAME, N, INT( CALLS(1) ), INT( CALLS(2) ),
     *                     STAT, F, CPU(1), CPU(2) 
      CLOSE( INPUT  )
      CALL CUTEST_uterminate( status )
      STOP

  910 CONTINUE
      WRITE( iout, "( ' CUTEst error, status = ', i0, ', stopping' )") 
     *   status
      STOP

  990 CONTINUE
      WRITE( out, "( ' Allocation error, status = ', I0 )" ) status
      STOP
C
C  Non-executable statements.
C
 1000 FORMAT( 18( D10.3, / ), 2( I10, / ), 7( L10, / ), L10 )
 1001 FORMAT( 16( D10.3, / ), 2( I10, / ), 7( L10, / ), L10 )
 2000 FORMAT( /, 24('*'), ' CUTEr statistics ', 24('*') //
     *    ,' Code used               :  CG_DESCENT',     /
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
C2020 FORMAT(  1X, A10, 1P, 2D12.4 )
 2030 FORMAT(  /, ' ** Warning from CGDMA. The problem as stated',
     *            ' includes simple bounds. ', /,
     *            '    These bounds will be ignored. ' )
C 2100 FORMAT( 1P, D10.3, 4X, 
C     * 'delta        Wolfe line search parameter',
C     * /, 1P, D10.3, 4X, 
C     * 'sigma        Wolfe line search parameter',
C     * /, 1P, D10.3, 4X, 
C     * 'epsilon      approximate Wolfe threshold factor',
C     * /, 1P, D10.3, 4X, 
C     * 'theta        update',
C     * /, 1P, D10.3, 4X, 
C     * 'gamma        required decay factor in interval',
C     * /, 1P, D10.3, 4X, 
C     * 'rho          growth factor in bracket',
C     * /, 1P, D10.3, 4X, 
C     * 'eta          lower bound for cg''s beta_k',
C     * /, 1P, D10.3, 4X, 
C     * 'psi0         factor used in very initial starting guess',
C     * /, 1P, D10.3, 4X, 
C     * 'psi1         factor previous step multiplied by in QuadStep',
C     * /, 1P, D10.3, 4X, 
C     * 'psi2         factor previous step is multipled by for startup')
C 2110 FORMAT( 1P, D10.3, 4X, 
C     * 'QuadCutOff   lower bound on rel change in f before QuadStep',
C     * /, 1P, D10.3, 4X, 
C     * 'restart_fac  restart cg in restart_fac*n iterations',
C     * /, 1P, D10.3, 4X, 
C     * 'maxit_fac    terminate in maxit_fac*n iterations',
C     * /, 1P, D10.3, 4X, 
C     * 'feps         stop when value change <= feps*|f|',
C     * /, I10, 4X, 
C     * 'nexpand      number of grow/shrink allowed in bracket',
C     * /, I10, 4X, 
C     * 'nsecant      number of secant steps allowed in line search',
C     * /, L10, 4X, 
C     * 'QuadStep     use initial quad interpolation in line search',
C     * /, L10, 4X, 
C     * 'PrintLevel   F (no print) T (intermediate results)')
C 2120 FORMAT(  L10, 4X, 
C     * 'PrintFinal   F (no print) T (print error messages, ',
C     * 'final error)',
C     * /, L10, 4X, 
C     * 'StopRule     F (|grad|_infty <= tol) T (... <= tol*(1+|f|))',
C     * /, L10, 4X, 
C     * 'ERule        F (eps_k = eps|f|) T (eps_k = eps)',
C     * /, L10, 4X, 
C     * 'AWolfe       F (Wolfe) T (+approx Wolfe) 2 (epsilon pert)',
C     * /, L10, 4X, 
C     * 'Step         F (no initial line search guess) T (guess in ',
C     * 'gnorm')
      END

      SUBROUTINE EVALF( F, X, N )

C  Evaluate the objective function

      INTEGER N
CS    REAL             F
CD    DOUBLE PRECISION F
CS    REAL             X( N )
CD    DOUBLE PRECISION X( N )
      EXTERNAL CUTEST_ufn

      INTEGER :: status

      CALL CUTEST_ufn( status, N, X, F )
      IF ( status /= 0 ) STOP

      RETURN
      END
      
      SUBROUTINE EVALG( G, X, N )

C  Evaluate the gradiuent of the objective function

      INTEGER N
CS    REAL             X( N ), G( N )
CD    DOUBLE PRECISION X( N ), G( N )
      EXTERNAL CUTEST_ugr

      INTEGER :: status

      CALL CUTEST_ugr( status, N, X, G )
      IF ( status /= 0 ) STOP

      RETURN
      END
