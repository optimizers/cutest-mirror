C     ( Last modified on 14 Jan 2013 at 15:00:00 )

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     Main program for SNOPTB using CUTEst
C
C     September 2004
C     CUTEst evolution, Nick Gould, January 2013
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      PROGRAM           SNOPT_main
      implicit none
      INTEGER :: m, n, ne, nb, nncon, nnjac, nnobj, nden, iobj, inform
      INTEGER :: ispecs, iprint, isumm, mincw, miniw, minrw, ns, ninf
      INTEGER :: lencw, leniu, leniw, lenru, lenrw, jslack, status
      INTEGER :: i, ii, indv, indf, j, jstrt, k, l_j, neq, nlc, njac
      INTEGER, PARAMETER :: input = 55, out = 6
      INTEGER, PARAMETER :: io_buffer = 11
      INTEGER, PARAMETER :: ldenj = 105
      DOUBLE PRECISION, PARAMETER :: zero = 0.0D+0, big = 1.0D+20
      DOUBLE PRECISION :: CPU( 2 ), CALLS( 7 )
      CHARACTER ( LEN = 8 ) :: start, probnm
      CHARACTER ( LEN = 10 ) :: pname
      INTEGER * 4, ALLOCATABLE, DIMENSION( : ) :: HA, HS
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: KA, IU, IW
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : ) :: X, BL, BU, RC
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : ) :: AA, Y, C, RU, RW
      LOGICAL, ALLOCATABLE, DIMENSION( : ) :: EQUATN, LINEAR
      CHARACTER ( LEN = 8 ), ALLOCATABLE, DIMENSION( : ) :: NAMES, CW
      CHARACTER ( LEN = 10 ), ALLOCATABLE, DIMENSION( : ) :: VNAME
      CHARACTER ( LEN = 10 ), ALLOCATABLE, DIMENSION( : ) :: CNAME
      DOUBLE PRECISION :: objadd, sinf, obj
      EXTERNAL :: SNOPT_evalcj, SNOPT_evalfg

C  Open the problem data file.

      OPEN ( input, FILE = 'OUTSDIF.d', FORM = 'FORMATTED',
     *       STATUS = 'OLD' )
      REWIND( input )

C  compute problem dimensions

      CALL CUTEST_cdimen( status, input, n, m )
      IF ( status /= 0 ) GO TO 910

C  allocate space 

      nb = n + m + 1
      lencw = 500
      leniw = 500
      lenrw = 500
      ALLOCATE( HS( nb ), KA( n + 1 ), X( nb ), BL( nb ), BU( nb ), 
     *          Y( m + 1 ), C( m + 1 ), RC( nb ), EQUATN( m + 1 ), 
     *          LINEAR( m + 1 ), VNAME( n ), CNAME( m + 1 ), 
     *          NAMES( nb ), CW( lencw ), IW( leniw ), RW( lenrw ), 
     *          STAT = status )
      IF ( status /= 0 ) GO TO 990

C  Set up the unit numbers for the SNOPT files.
C  ispecs is the Specifications file.
C  iprint is the Print file.
C  isumm is the Summary file.
C
      ispecs = 4
      iprint = 15
      isumm  = 6
      CALL S1OPEN( ispecs, 1, 'IN ' )
      CALL S1OPEN( iprint, 2, 'OUT' )

C     Set options to default values and read Specs file.
C
C     ------------------------------------------------------------------
C     Set all options as undefined.
C     ------------------------------------------------------------------

      CALL SNINIT( iprint, isumm, CW, lencw, IW, leniw, RW, lenrw )

C     ------------------------------------------------------------------
C     Read a Specs file.
C     ------------------------------------------------------------------

      CALL SNSPEC( ispecs, inform, CW, lencw, IW, leniw, RW, lenrw )
      IF ( inform .GE. 2 ) THEN
         IF ( out .GT. 0 ) WRITE( out, 2010 )
         STOP
      END IF
      nden = IW( ldenj )
      nden = MAX( nden, 1 )

C  input problem data using csetup

      CALL CUTEST_csetup( status, input, out, io_buffer, n, m, 
     *                    X, BL, BU, Y, BL( n + 1 ), BU( n + 1 ), 
     *                    EQUATN, LINEAR, 0, 1, 1 )
      CLOSE( input )
      IF ( status /= 0 ) GO TO 910

C  compute the numbers of nonlinear variables, and linear/equatity constraints

      CALL CUTEST_cstats( status, nnobj, nnjac, neq, nlc )
      IF ( status /= 0 ) GO TO 910

!  compute the objective and constraints at X = 0

      DO 90 i = 1, n
        RC( i ) = zero
   90 CONTINUE
      CALL CUTEST_cfn( status, n, m, RC, obj, C )
      IF ( status /= 0 ) GO TO 910

C  Determine the number of nonlinear constraints

      nncon = m - nlc

C  Use the constraint bounds to set the bounds on the slack variables

      DO 100 i = 1, m
        IF ( EQUATN( i ) ) THEN
          BL( n + i ) = zero
          BU( n + i ) = zero
        END IF
  100 CONTINUE

C  Add one to m for linear objective row. Also set BL and BU for the objective 
C  row. If the objective function has a linear part, set iobj to m

      m = m + 1
      BL( n + m ) = - big
      BU( n + m ) = big
      C( m ) = zero
C     X( n + m ) = zero
      IF ( nnobj .LT. n ) THEN
        iobj = m
      ELSE
        iobj = 0
      END IF

C  Set up AA(i), KA(j) and HA(i). AA(i) gives the i-th element in the Jacobian.
C  KA(j) gives starting address in AA of entries for variable j. HA(i) gives 
C  constraint index for i-th element in AA

C  Jacobian is to be stored in dense format

      IF ( nden .EQ. 1 ) THEN

C  Allocate space for AA and HA as well as user workspace

        leniu = 1
        lenru = n
        ne = m * n
        ALLOCATE( HA( ne ), AA( ne ), IU( leniu ), RU( lenru ), 
     *            STAT = status )
        IF ( status /= 0 ) GO TO 990

        IU( 1 ) = nden

C  find the entries in the dense Jacobian

        CALL CUTEST_cgr( status, n, m, X, Y, .FALSE., 
     *                   RU, .FALSE., m, n, AA )
        IF ( status /= 0 ) GO TO 910

C  Set KA(j) and HA(j)

        IF ( ne .GT. 0 ) THEN
          KA( 1 ) = 1
          DO 220 j = 1, n 
            KA( j + 1 ) = KA( j ) + m
            k = KA( j ) - 1
            DO 210 i = 1, m
              HA( k + i ) = i
  210       CONTINUE

C  Copy gradient of linear part of objective function into row iobj of Jacobian

            IF ( j .GT. nnobj ) AA( k + iobj ) = RU( j )
  220     CONTINUE
        END IF

C  Jacobian is to be stored in sparse format

      ELSE

C  compute the number of nonzeros in the Jacobian

        CALL CUTEST_cdimsj( status, ne )
        IF ( status /= 0 ) GO TO 910

C  Partition the integer sparse work vector IU and RU

        jstrt = 4
        indv = jstrt + n + 1
        indf = indv + ne
        leniu = indf + ne
        lenru = ne
        l_j = ne

C  Allocate space for AA and HA as well as user workspace

        ALLOCATE( HA( ne ), AA( ne ), IU( leniu ), RU( lenru ), 
     *            STAT = status )
        IF ( status /= 0 ) GO TO 990

        IU( 1 ) = nden
        IU( 2 ) = jstrt
        IU( 3 ) = indv
        IU( 4 ) = indf

C  Use CSGR to find entries in sparse Jacobian. Since CSGR and MINOS use 
C  different sparse formats, store Jacobian temporarily in RU and IU.

        CALL CUTEST_csgr( status, n, m, X, Y, .FALSE., ne, l_j,
     *                    RU, IU( indv + 1 ), IU( indf + 1 ) )
        IF ( status /= 0 ) GO TO 910
        k = ne

C  Initialize KA

        DO 250 j = 1, n
          KA( j ) = 0
  250   CONTINUE

C  Count Jacobian entries for each variable j. Store counts in KA(j).
C  Don't include nonlinear objective function entries

        DO 300 ii = 1, ne
          j = IU( indv + ii )
          i = IU( indf + ii )
          IF ( i .GT. 0 .OR. j .GT. nnobj ) THEN
            KA( j ) = KA( j ) + 1
          ELSE
            k = k - 1
          END IF
  300   CONTINUE
        KA( n + 1 ) = k + 1 

C  Now set KA(j) to starting address for variable j

        DO 310 j = n, 1, - 1
          KA ( j ) = KA( j + 1 ) - KA( j )
          IU( jstrt + j ) = 0
  310   CONTINUE

C  Loop through nonlinear Jacobian entries. Put correct entries in AA and HA.
C  Use KA to keep track of position for each variable j. Also count nonlinear 
C  Jacobian entries for each variable j. Store count in IU(jstrt+j)

        njac = 0
        DO 320 k = 1, ne
          j = IU( indv + k )
          i = IU( indf + k )
          IF ( i .GT. 0 .AND. i .LE. nncon .AND. j .LE. nnjac ) THEN
            ii = KA( j )
            AA( ii ) = RU( k )
            HA( ii ) = i
            KA( j ) = ii + 1
            IU( jstrt + j ) = IU( jstrt + j ) + 1
            njac = njac + 1
          END IF
  320   CONTINUE

C  Now loop through linear Jacobian entries, including linear objective 
C  function entries. Put correct entries in AA and HA. Use KA to keep track 
C  of position for each variable J.

        DO 330 k = 1, ne
          j = IU( indv + k )
          i = IU( indf + k )
          IF ( i .EQ. 0 .AND. j .GT. nnobj ) THEN
            ii = KA ( j )
            AA( ii ) = RU( k )
            HA( ii ) = iobj
            KA( j ) = ii + 1
          ELSE IF ( ( i .GT. 0 .AND. j .GT. nnjac ) .OR.
     *                i .GT. nncon ) THEN
            ii = KA ( j )
            AA( ii ) = RU( k )
            HA( ii ) = i
            KA( j ) = ii + 1
          END IF
  330   CONTINUE

C  Reset KA(j) and set IU(jstrt+j). IU(jstrt+j) now gives starting address for
C  variable j in nonlinear Jacobian. These addresses are needed in SNOPT_evalcj

        IU( jstrt + n + 1 ) = njac + 1
        DO 340 j = n, 2, - 1
          KA ( j ) = KA( j - 1 )
          IU( jstrt + j ) = IU( jstrt + j + 1 ) - IU( jstrt + j )
  340   CONTINUE
        KA( 1 ) = 1 
        IU( jstrt + 1 ) = 1
        ne = KA( n + 1 ) - 1
      END IF

C  Incorporate nonzero constants from linear constraints as bounds on slack 
C  variables - constants for nonlinear constraints are added in CUTEST_ccfg 
C  or  CUTEST_ccfsg, which are called by SNOPT_evalcj

      DO 400 i = 1, m
        jslack = n + i
        IF ( i .GT. nncon ) THEN
          BU( jslack ) = BU( jslack ) + C( i )
          BL( jslack ) = BL( jslack ) + C( i )
        END IF

C  If possible, set slack variables to be nonbasic at zero.

        X( jslack ) = MAX( zero, BL( jslack ) )
        X( jslack ) = MIN( X( jslack ), BU( jslack ) )
  400 CONTINUE

C  Incorporate nonzero constants from objective function groups only if 
C  objective function is completely linear - iIf the objective function has 
C  nonlinear part, constants are added in COFG, which is called by SNOPT_evalfg

      objadd = zero
      IF ( nnobj .EQ. 0 ) objadd = objadd - obj

C  determine the names for problem quantities

      CALL CUTEST_cnames( status, n, m, pname, VNAME, CNAME )
      IF ( status /= 0 ) GO TO 910

C  Assign SNOPT names to the problem, variables and constraints

      probnm = pname( 1 : 8 )
      DO 410 j = 1, n
        NAMES( j ) = VNAME( j )( 1 : 8 )
  410 CONTINUE
      DO 420 i = 1, m - 1
        NAMES( n + i ) = CNAME( i )( 1 : 8 )
  420 CONTINUE
      NAMES( n + m ) = probnm

      nb = n + m
      DO 510 I = 1, nb
        HS( i ) = 0
  510 CONTINUE
      DO 520 I = 1, nncon
        Y( i ) = zero
  520 CONTINUE

C  Call SNOPT as a subroutine

  590 CONTINUE
      start = 'COLD'
      CALL SNOPTB( start, m, n, ne, nb, nncon, nnobj, nnjac, iobj, 
     *             objadd, probnm, SNOPT_evalcj, SNOPT_evalfg,
     *             AA, HA, KA, BL, BU, NAMES, HS, X, Y, RC, inform,
     *             mincw , miniw, minrw, ns, ninf, sinf, obj,
     *             CW, lencw, IU, leniu, RU, lenru, 
     *             CW, lencw, IW, leniw, RW, lenrw )

C  if there is insufficient workspace, try increasing it

      IF ( inform .GE. 81 .AND. inform .LE. 84 ) THEN
        DEALLOCATE( CW, IW, RW, STAT = status )
        IF ( status /= 0 ) GO TO 990
        lencw = MAX( lencw, 3 * mincw / 2 )
        leniw = MAX( leniw, 3 * miniw / 2 )
        lenrw = MAX( lenrw, 3 * minrw / 2 )
        ALLOCATE( CW( lencw ), IW( leniw ), RW( lenrw ), STAT = status )
        IF ( status /= 0 ) GO TO 990
        GO TO 590
      END IF

C  Try to handle abnormal SNOPT inform codes gracefully

      IF ( inform .GE. 20 .AND.
     *  ( iprint .GT. 0 .OR. isumm .GT. 0 ) ) THEN
         IF ( iprint .GT. 0 ) WRITE ( iprint, 3000 ) inform
         IF ( isumm .GT. 0 ) WRITE ( isumm, 3000 ) inform
         IF ( inform .EQ. 42  .OR. inform .EQ. 43 ) THEN
            IF ( iprint .GT. 0 ) WRITE( iprint, 3010 )
            IF ( isumm .GT. 0 ) WRITE( isumm , 3010 )
         END IF
      END IF
      CALL CUTEST_creport( status, CALLS, CPU )
      IF ( status /= 0 ) GO TO 910
      WRITE ( out, 2020 ) NAMES( 1 ), n, m - 1, CALLS( 1 ), CALLS( 2 ), 
     *    CALLS( 5 ), CALLS( 6 ), inform, obj, CPU( 1 ), CPU( 2 )
      CLOSE( ispecs )
      CLOSE( iprint )
      DEALLOCATE( HS, KA, X, BL, BU, Y, C, RC, EQUATN, LINEAR, 
     *            VNAME, CNAME, NAMES, HA, AA, CW, IU, IW, RU, RW, 
     *            STAT = status )
      CALL CUTEST_cterminate( status )
      STOP

  910 CONTINUE
      WRITE( out, "( ' CUTEst error, status = ', i0, ', stopping' )") 
     *   status
      STOP

  990 CONTINUE
      WRITE( out, "( ' Allocation error, status = ', I0 )" ) status
      STOP

C  Non-executable statements

C2000 FORMAT( /, ' ** SUBROUTINE SNOPT_main: array length ', A, 
C    *        ' too small.', /, ' -- Minimization abandoned.',
C    *        /, ' -- Increase the parameter ', A, ' by at least ', I0,
C    *           ' and restart.'  )
 2010 FORMAT( /, ' ** PROGRAM SNOPT_main: No Specs file found.' )
 2020 FORMAT( /, 24('*'), ' CUTEst statistics ', 24('*') //
     *    ,' Code used               :  SNOPT',    /
     *    ,' Problem                 :  ', A10,    /
     *    ,' # variables             =      ', I10 /
     *    ,' # constraints           =      ', I10 /
     *    ,' # objective functions   =        ', F8.2 /
     *    ,' # objective gradients   =        ', F8.2 / 
     *    ,' # constraints functions =        ', F8.2 /
     *    ,' # constraints gradients =        ', F8.2 /
     *     ' Exit code               =      ', I10 /
     *    ,' Final f                 = ', E15.7 /
     *    ,' Set up time             =      ', 0P, F10.2, ' seconds' /
     *     ' Solve time              =      ', 0P, F10.2, ' seconds' //
     *     66('*') / )
 3000 FORMAT( /, ' WARNING!  Abnormal SNOPT termination code:',
     *           ' inform = ', I2 )
 3010 FORMAT(    ' Not enough storage to solve the problem.',
     *        /, ' Reduce parameters in SNOPT.SPC file or increase',
     *           ' LENIW and LENRW in SNOPT_main.' )
      END

C*********************************************************************
C
C SUBROUTINE  SNOPT_evalfg
C
C*********************************************************************

      SUBROUTINE SNOPT_evalfg( mode, n, X, f, G, nstate, CW, lencw, 
     *                         IW, leniw, RW, lenrw )
      INTEGER :: mode, n, nstate, lencw, leniw, lenrw
      INTEGER :: IW( leniw )
      DOUBLE PRECISION :: f
      DOUBLE PRECISION  :: X( n ), G( n ), RW( lenrw )
      CHARACTER ( LEN = 8 ) :: CW( lencw )

C  Local variables

      INTEGER :: status
      LOGICAL :: grad
      IF ( mode .EQ. 0 ) THEN
        grad = .FALSE.
      ELSE
        grad = .TRUE.
      END IF
      CALL CUTEST_cofg( status, n, X, f, G, grad )
      IF ( status .NE. 0 ) THEN
        WRITE( 6, "( ' CUTEst error, status = ', i0, ', stopping' )") 
     *   status
        STOP
      END IF
      RETURN
      END

C*********************************************************************
C
C SUBROUTINE  SNOPT_evalcj
C
C*********************************************************************

      SUBROUTINE SNOPT_evalcj( mode, m, n, njac, X, C, JAC, nstate,
     *                         CU, lencu, IU, leniu, RU, lenru )
      INTEGER :: mode, m, n, njac, nstate, lencu, leniu, lenru
      INTEGER :: IU( leniu )
      DOUBLE PRECISION :: X( n ), C( m ), JAC( njac ), RU( lenru )
      CHARACTER ( LEN = 8 ) :: CU( lencu )

C  Local variables

      INTEGER :: i, j, k, indv, indf, jstrt, nnzj, status
      LOGICAL :: grad

      IF ( mode .EQ. 0 ) THEN
        grad = .FALSE.
      ELSE
        grad = .TRUE.
      END IF

C  Jacobian is stored in dense format

      IF ( IU( 1 ) .EQ. 1 ) THEN
        CALL CUTEST_ccfg( status, n, m, X, C, .FALSE., m, n, JAC, grad )
        IF ( status .NE. 0 ) GO TO 910

C  Jacobian is stored in sparse format

      ELSE
        jstrt = IU( 2 )
        indv = IU( 3 )
        indf = IU( 4 )
        CALL CUTEST_ccfsg( status, n, m, X, C, nnzj, njac, RU, 
     *                     IU( indv + 1 ), IU( indf + 1 ), grad )
        IF ( status .NE. 0 ) GO TO 910

C  Copy Jacobian from CCFSG, contained in RU, into MINOS Jacobian JAC in 
C  correct order. Use IU(jstrt+j) to keep track of position for variable j

        IF ( grad ) THEN
          DO 130 i = 1, nnzj
            j = IU( indv + i )
            k = IU( jstrt + j )
            JAC( k ) = RU( i )
            IU( jstrt + j ) = K + 1
  130     CONTINUE

C  Reset IU(jstrt+j)

          DO 140 j = n, 2, - 1
            IU( jstrt + j ) = IU( jstrt + j - 1 )
  140     CONTINUE
          IU( jstrt + 1 ) = 1
        END IF
      END IF
      RETURN

  910 CONTINUE
      WRITE( 6, "( ' CUTEst error, status = ', i0, ', stopping' )") 
     *   status
        STOP
      END
