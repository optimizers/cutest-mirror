C     ( Last modified on 11 Jan 2013 at 14:30:00 )

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     Main program for MINOS using CUTEst
C
C     Ingrid Bongartz, August 1992
C     CUTEst evolution, Nick Gould, January 2013
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      PROGRAM           MINOS_main
      INTEGER           NWCORE
CBIG  PARAMETER       ( NWCORE = 1000000 )
CMED  PARAMETER       ( NWCORE = 200000 )
CTOY  PARAMETER       ( NWCORE = 50000 )
CCUS  PARAMETER       ( NWCORE = 500000 )
      DOUBLE PRECISION  Z( NWCORE )
      INTEGER :: ispecs, iprint, isumm, ns, nin, i, nnden
      INTEGER :: m, n, ne, nb, nncon, nnjac, nnobj, iobj, inform, mincor
      INTEGER, PARAMETER :: input = 55, out = 6
      DOUBLE PRECISION :: objadd, sinf, obj
      DOUBLE PRECISION, PARAMETER :: zero = 0.0D+0
      DOUBLE PRECISION :: CPU( 2 ), CALLS( 7 )

      INTEGER           MMAX, NMAX, NEMAX, NBMAX 
CBIG  PARAMETER       ( MMAX = 5000 )
CMED  PARAMETER       ( MMAX = 500 )
CTOY  PARAMETER       ( MMAX = 100 )
CCUS  PARAMETER       ( MMAX = 1000 )
      PARAMETER       ( NMAX = 3 * MMAX, NEMAX = 5 * NMAX,
     *                  NBMAX = NMAX + MMAX )
      INTEGER           NAME1( 1 ), NAME2 ( 1 )
      DOUBLE PRECISION  AA   ( NEMAX ), X     ( NBMAX ), BL( NBMAX ),
     *                  BU   ( NBMAX ), V     ( MMAX  ), RC( NBMAX )
      INTEGER * 4       HA   ( NEMAX ), HS    ( NBMAX )
      INTEGER           KA   ( NMAX + 1 )
      LOGICAL           EQUATN( MMAX ), LINEAR( MMAX )
      CHARACTER * 10    VNAME( NMAX ), CNAME( MMAX )
      CHARACTER * 8     START, NAMES( 5 )

C  MINOS common block

      INTEGER ::  ncom, nden, nlag, nmajor, nminor
      DOUBLE PRECISION  :: penpar, rowtol
      COMMON   /M8AL1 / penpar, rowtol, ncom, nden, nlag, nmajor, nminor

C  Open the relevant file.

      OPEN ( input, FILE = 'OUTSDIF.d', FORM = 'FORMATTED',
     *       STATUS = 'OLD' )
      REWIND( input )

C  Set up the unit numbers for the MINOS files:
C  ispecs is the Specifications file
C  iprint is the Print file
C  isumm is the Summary file

      ispecs = 4
      iprint = 9
      isumm  = 0
      CALL m1open( ispecs, 1, 'IN' )
      CALL m1open( iprint, 2, 'OUT' )

C  Set options to default values and read Specs file

      CALL MISPEC( ispecs, iprint, isumm, nwcore, inform )
      IF ( inform .GE. 2 ) THEN
        IF ( out .GT. 0 ) WRITE( out, 2010 )
        STOP
      END IF
      nnden = nden

C  input problem data using MINOS_setup, which calls CSETUP

      CALL MINOS_setup( input, out, MMAX, NMAX, NEMAX,
     *             NBMAX, NNDEN, M, N, NE,
     *             NNCON, NNJAC, NNOBJ, IOBJ,
     *             OBJADD, HA, KA, AA,
     *             X, BL, BU, V,
     *             NAMES, VNAME, CNAME, EQUATN, LINEAR )
      NB = N + M
      DO 10 I = 1, NB
         HS( I ) = 0
         RC( I ) = zero
   10 CONTINUE

C  Call MINOS as a subroutine

      START = 'COLD'
      CALL MINOSS( START, M, N, NB, NE, 1,
     *             NNCON, NNOBJ, NNJAC,
     *             IOBJ, OBJADD, NAMES,
     *             AA, HA, KA, BL, BU,
     *             NAME1, NAME2, 
     *             HS, X, V, RC,
     *             INFORM, MINCOR, NS, NINF, SINF,
     *             OBJ, Z, NWCORE )
      CALL CUTEST_creport( status, CALLS, CPU )

C  Try to handle abnormal MINOS inform codes gracefully

      IF ( INFORM .GE. 20 .AND.
     *  ( IPRINT .GT. 0 .OR. ISUMM .GT. 0 ) ) THEN
         IF ( IPRINT .GT. 0 ) WRITE ( IPRINT, 3000 ) INFORM
         IF ( ISUMM  .GT. 0 ) WRITE ( ISUMM, 3000 ) INFORM
         IF ( INFORM .EQ. 20 ) THEN
            IF ( IPRINT .GT. 0 ) WRITE ( IPRINT, 3020 )
            IF ( ISUMM  .GT. 0 ) WRITE ( ISUMM, 3020 )
         ELSE IF ( INFORM .EQ. 21 ) THEN
            IF ( IPRINT .GT. 0 ) WRITE ( IPRINT, 3021 )
            IF ( ISUMM  .GT. 0 ) WRITE ( ISUMM, 3021 )
         ELSE IF ( INFORM .EQ. 22 ) THEN
            IF ( IPRINT .GT. 0 ) WRITE ( IPRINT, 3022 )
            IF ( ISUMM  .GT. 0 ) WRITE ( ISUMM, 3022 )
         ELSE IF ( INFORM .EQ. 32 ) THEN
            IF ( IPRINT .GT. 0 ) WRITE ( IPRINT, 3032 )
            IF ( ISUMM  .GT. 0 ) WRITE ( ISUMM, 3032 )
         ELSE IF ( INFORM .EQ. 42 ) THEN
            IF ( IPRINT .GT. 0 ) WRITE ( IPRINT, 3042 )
            IF ( ISUMM  .GT. 0 ) WRITE ( ISUMM, 3042 )
         END IF
      END IF
      IF ( out .GT. 0 ) THEN
         WRITE( out, 2110 ) OBJ, ( I, VNAME( I ), X( I ), BL( I ), 
     *                    BU( I ), RC( I ), I = 1, N )
         IF ( M .GT. 1 ) WRITE( out, 2120 ) ( I, CNAME( I ), 
     *        X( N + I ), BL( N + I ), BU( N + I ), V( I ), 
     *        EQUATN( I ), LINEAR( I ), I = 1, M - 1 )
      END IF
      IF ( ISUMM  .GT. 0 ) 
     *  WRITE ( ISUMM, 2020 ) NAMES(1), N, M, CALLS(1), 
     *                         CALLS(2), CALLS(5), CALLS(6), INFORM,
     *                         OBJ, CPU(1), CPU(2)
      CLOSE( ISPECS )
      CLOSE( IPRINT )
      STOP

C  Non-executable statements

 2010 FORMAT( /, ' ** PROGRAM MINOSMA: No Specs file found.' )
 2020 FORMAT( /, 24('*'), ' CUTEst statistics ', 24('*') //
     *,' Code used               :  MINOS',    /
     *,' Problem                 :  ', A10,    /
     *,' # variables             =    ', I10 /
     *,' # constraints           =    ', I10 /
     *,' # objective functions   =      ', F8.2 /
     *,' # objective gradients   =      ', F8.2 / 
     *,' # constraints functions =      ', F8.2 /
     *,' # constraints gradients =      ', F8.2 /
     *,' Exit code               =    ', I10 /
     *,' Final f                 = ', E15.7 /
     *,' Set up time             =    ', 0P, F10.2, ' seconds' /
     *,' Solve time              =    ', 0P, F10.2, ' seconds' //
     *     66('*') / )
 2110 FORMAT( /, ' the objective function value: ', 1P, D12.4, /,
     *        /, ' the variables:', //,
     *        '     # name          value    lower bound upper bound', 
     *        '  dual value', 
     *        /, ( I6, 1X, A10, 1P, 4D12.4 ) )
 2120 FORMAT( /, ' the constraints:', //,
     *        '     # name          value    lower bound upper bound', 
     *        '  multiplier equality linear', 
     *        /, ( I6, 1X, A10, 1P, 4D12.4, 4X, L1, 8X, L1 ) )
 3000 FORMAT( /, ' WARNING!  Abnormal MINOS termination code:',
     *           ' INFORM = ', I2 )
 3020 FORMAT(    ' Not enough storage for the basis factorization.',
     *        /, ' Reduce parameters in MINOS.SPC file or increase',
     *           ' NWCORE in MINOSMA.' )
 3021 FORMAT(    ' Error in basis package.' )
 3022 FORMAT(    ' The basis is singular after several attempts to',
     *        /, ' factorize it (and add slacks where necessary).' )
 3032 FORMAT(    ' System error.  Wrong number of basic variables.' )
 3042 FORMAT(    ' Not enough storage to solve the problem.',
     *        /, ' Reduce parameters in MINOS.SPC file or increase',
     *           ' NWCORE in MINOSMA.' )
      END

      SUBROUTINE MINOS_setup( input, out, MMAX, NMAX, NEMAX,
     *                   NBMAX, NDEN, M, N, NE,
     *                   NNCON, NNJAC, NNOBJ, IOBJ,
     *                   OBJADD, HA, KA, AA,
     *                   X, BL, BU, V,
     *                   NAMES, VNAME, CNAME, EQUATN, LINEAR )
      INTEGER :: input, out, mmax, nmax, nemax, nbmax, nden, m, n
      INTEGER :: ne, nncon, nnjac, nnobj, iobj
      INTEGER :: HA( nemax ), KA( nmax + 1 )
      DOUBLE PRECISION :: objadd
      DOUBLE PRECISION :: AA( nemax ), X ( nbmax ), BL( nbmax ), 
     *                    BU( nbmax ), V ( mmax )
      CHARACTER ( LEN = 8 ) :: NAMES ( 5 )
      CHARACTER ( LEN = 10 ) :: pname, VNAME( nmax ), CNAME( mmax )
      LOGICAL :: EQUATN( mmax  ), LINEAR( mmax  )
C
C  Set up the input data for MINOS.
C
C  Modified November 1992 to call CSETUP.
C
      INTEGER            LIWK, LWK, LFUVAL, LLOGIC, LCHARA
C
C  ---------------------------------------------------------------------
C
C  Parameters whose value might be changed by the user:
C
C  The following parameters define the sizes of problem
C  dependent arrays. These may be changed by the user to
C  suit a particular problem or system configuration.
C
C  The TOOLS will issue error messages if any of these sizes
C  is too small, telling which parameter to increase.
C
C  ---------------------------------------------------------------------
C
C#{sizing}
      INTEGER                IWK( LIWK    )
      DOUBLE PRECISION       WK ( LWK     )
      LOGICAL              LOGI ( LLOGIC  )
      CHARACTER * 10       CHA  ( LCHARA  )
      DOUBLE PRECISION  FUVALS  ( LFUVAL  )
C
C  ---------------------------------------------------------------------
C
C  End of parameters which might be changed by the user.
C
C  ---------------------------------------------------------------------
C
C  integer variables from the GLOBAL common block.
C
      INTEGER            NG, NELNUM, NGEL, NVARS, NNZA, NGPVLU
      INTEGER            NEPVLU, NG1, NEL1, ISTADG, ISTGP, ISTADA
      INTEGER            ISTAEV, ISTEP, ITYPEG, KNDOFC, ITYPEE
      INTEGER            IELING, IELVAR, ICNA, ISTADH, INTVAR, IVAR
      INTEGER            ICALCF, ITYPEV, IWRK, A, B     
      INTEGER            U, GPVALU, EPVALU
      INTEGER            ESCALE, GSCALE, VSCALE, GVALS, XT, DGRAD 
      INTEGER            Q, WRK, INTREP, GXEQX, GNAMES, VNAMES
      INTEGER            LO, CH, LIWORK, LWORK, NGNG, FT    
      INTEGER            LA, LB, NOBJGR, LU, LELVAR
      INTEGER            LSTAEV, LSTADH, LNTVAR, LCALCF
      INTEGER            LELING, LINTRE, LFT, LGXEQX, LSTADG, LGVALS
      INTEGER            LICNA, LSTADA, LKNDOF, LGPVLU, LEPVLU
      INTEGER            LGSCAL, LESCAL, LVSCAL, LCALCG
C
C  integer variables from the LOCAL common block.
C
      INTEGER            LFXI, LGXI, LHXI, LGGFX, LDX, LGRJAC
      INTEGER            LQGRAD, LBREAK, LP, LXCP, LX0, LGX0  
      INTEGER            LDELTX, LBND, LWKSTR, LSPTRS, LSELTS, LINDEX
      INTEGER            LSWKSP, LSTAGV, LSTAJC, LIUSED, LFREEC
      INTEGER            LNNONZ, LNONZ2, LSYMMD, LSYMMH
      INTEGER            LSLGRP, LSVGRP, LGCOLJ, LVALJR, LSEND 
      INTEGER            LNPTRS, LNELTS, LNNDEX, LNWKSP, LNSTGV
      INTEGER            LNSTJC, LNIUSE, LNFREC, LNNNON, LNNNO2, LNSYMD
      INTEGER            LNSYMH, LNLGRP, LNVGRP, LNGCLJ, LNVLJR, LNQGRD
      INTEGER            LNBRAK, LNP, LNBND, LNFXI, LNGXI, LNGUVL
      INTEGER            LNHXI, LNHUVL, LNGGFX, LNDX, LNGRJC, LIWK2 
      INTEGER            LWK2, MAXSIN, NINVAR, MAXSEL
      INTEGER            NTYPE, NSETS, LSTYPE, LSSWTR, LSSIWT, LSIWTR
      INTEGER            LSWTRA, LNTYPE, LNSWTR, LNSIWT, LNIWTR
      INTEGER            LNWTRA, LSISET, LSSVSE, LNISET, LNSVSE
      LOGICAL            ALTRIV, FIRSTG
      COMMON / GLOBAL /  IWK, WK, FUVALS, LOGI,  
     *                   NG, NELNUM, NGEL, NVARS, NNZA, NGPVLU,
     *                   NEPVLU, NG1, NEL1, ISTADG, ISTGP, ISTADA,
     *                   ISTAEV, ISTEP, ITYPEG, KNDOFC, ITYPEE, 
     *                   IELING, IELVAR, ICNA, ISTADH, INTVAR, IVAR,
     *                   ICALCF, ITYPEV, IWRK, A, B, 
     *                   U, GPVALU, EPVALU, 
     *                   ESCALE, GSCALE, VSCALE, GVALS, XT, DGRAD,
     *                   Q, WRK, INTREP, GXEQX, GNAMES, VNAMES,
     *                   LO, CH, LIWORK, LWORK, NGNG, FT,
     *                   ALTRIV, FIRSTG,
     *                   LA, LB, NOBJGR, LU, LELVAR,
     *                   LSTAEV, LSTADH, LNTVAR, LCALCF,
     *                   LELING, LINTRE, LFT, LGXEQX, LSTADG, LGVALS,
     *                   LICNA, LSTADA, LKNDOF, LGPVLU, LEPVLU,
     *                   LGSCAL, LESCAL, LVSCAL, LCALCG
      COMMON / CHARA /   CHA
      COMMON / LOCAL /   LFXI, LGXI, LHXI, LGGFX, LDX, LGRJAC,
     *                   LQGRAD, LBREAK, LP, LXCP, LX0, LGX0,
     *                   LDELTX, LBND, LWKSTR, LSPTRS, LSELTS, LINDEX,
     *                   LSWKSP, LSTAGV, LSTAJC, LIUSED, LFREEC,
     *                   LNNONZ, LNONZ2, LSYMMD, LSYMMH,
     *                   LSLGRP, LSVGRP, LGCOLJ, LVALJR, LSEND,
     *                   LNPTRS, LNELTS, LNNDEX, LNWKSP, LNSTGV,
     *                   LNSTJC, LNIUSE, LNFREC, LNNNON, LNNNO2, LNSYMD,
     *                   LNSYMH, LNLGRP, LNVGRP, LNGCLJ, LNVLJR, LNQGRD,
     *                   LNBRAK, LNP, LNBND, LNFXI, LNGXI, LNGUVL,
     *                   LNHXI, LNHUVL, LNGGFX, LNDX, LNGRJC, LIWK2,
     *                   LWK2, MAXSIN, NINVAR, MAXSEL, NTYPE,
     *                   NSETS, LSTYPE, LSSWTR, LSSIWT, LSIWTR,
     *                   LSWTRA, LNTYPE, LNSWTR, LNSIWT, LNIWTR, 
     *                   LNWTRA, LSISET, LSSVSE, LNISET, LNSVSE
      INTEGER            NNOV, NNJV 
      COMMON / NNVARS /  NNOV, NNJV 
      SAVE             / GLOBAL /, / LOCAL /, / CHARA /, / NNVARS /
C
C  Sparse Jacobian common block 
C
      INTEGER            LICWK, LCWK
C    to accommodate increased MMAX for large installations.
CBIG  PARAMETER        ( LICWK = 165001 )
CMED  PARAMETER        ( LICWK = 16501 )
CTOY  PARAMETER        ( LICWK = 3301 )
CCUS  PARAMETER        ( LICWK = 165001 )
      INTEGER            ICWK( LICWK )
C    to accommodate increased MMAX for large installations.
CBIG  PARAMETER        ( LCWK  = 75000 )
CMED  PARAMETER        ( LCWK  = 7500 )
CTOY  PARAMETER        ( LCWK  = 1500 )
CCUS  PARAMETER        ( LCWK  = 75000 )
CD    DOUBLE PRECISION   CWK ( LCWK )
      INTEGER            JSTRT, INDV, INDF 
      COMMON / SPJAC /   CWK, ICWK, JSTRT, INDV, INDF
      SAVE             / SPJAC /
C
C  Local variable declarations
C
      INTEGER            I, II, IG, J, JG, K, JSLACK, MEND, NJAC
      LOGICAL            LTEMP
      DOUBLE PRECISION   ATEMP
      DOUBLE PRECISION, PARAMETER :: zero = 0.0D+0, big = 1.0D+20
C
C  input problem data using csetup.
C

      CALL CUTEST_csetup( status, input, out, N, M, X, BL, BU, 
     *             V, BL( NMAX + 1 ), BU( NMAX + 1 ), EQUATN, LINEAR, 
     *             0, 0, 1 )
      CLOSE( input )
      CALL CUTEST_cnames( status, N, M, PNAME, VNAME, CNAME )
      NNOBJ = NNOV
      NNJAC = NNJV
C
C  Determine the number of nonlinear constraints.  
C  Use the constraint bounds to set the bounds
C  on the slack variables.
C
      NNCON = 0
      DO 100 IG = 1, NG
         I = IWK( KNDOFC + IG )
         IF ( I .GT. 0 ) THEN
            IWK( LSEND + I ) = IG
            IF ( .NOT. LINEAR( I ) ) NNCON = NNCON + 1
            IF ( EQUATN( I ) ) THEN
               BL( N + I ) = zero
               BU( N + I ) = zero
            ELSE
               IF ( NMAX .GT. N ) THEN 
                  BL( N + I ) = - BU( NMAX + I )
                  BU( N + I ) = - BL( NMAX + I )
               ELSE
                  ATEMP       = - BU( N + I )
                  BU( N + I ) = - BL( N + I )
                  BL( N + I ) = ATEMP
               END IF
            END IF
         END IF
  100 CONTINUE
      IF ( NNCON .EQ. 0 .OR. NNCON .EQ. M ) GO TO 130
C
C  Reorder the constraints so that the nonlinear constraints occur before 
C  the linear ones.
C
      MEND = M
C
C  Run forward through the constraints until a linear constraint 
C  is encountered.
C
      DO 120 I = 1, M
         IF ( I .GT. MEND ) GO TO 130
         IG = IWK( LSEND + I )
C           write(6,*) ' group ', IG, ' type ', I, ' linear? ',
C     *                  LINEAR( I )
         IF ( LINEAR( I ) ) THEN
C
C  Constraint I is linear. Now, run backwards through the 
C  constraints until a nonlinear one is encountered.
C
            DO 110 J = MEND, I, - 1
               JG    = IWK( LSEND + J )
C              write(6,*) ' group ', JG, ' type ', J,
C     *                   ' linear? ', LINEAR( J )
               IF ( .NOT. LINEAR( J ) ) THEN
C                 write(6,*) ' swapping constraints ', I,
C     *                      ' and ', J
                  MEND = J - 1
C
C  Interchange the data for constraints I and J.
C
                  IWK ( LSEND  +  I ) = JG 
                  IWK ( LSEND  +  J ) = IG
                  IWK ( KNDOFC + IG ) = J
                  IWK ( KNDOFC + JG ) = I
                  LTEMP               = LINEAR( I )
                  LINEAR( I )         = LINEAR( J )
                  LINEAR( J )         = LTEMP
                  LTEMP               = EQUATN( I )
                  EQUATN( I )         = EQUATN( J )
                  EQUATN( J )         = LTEMP

C    estimates
                  ATEMP               = V     ( I )
                  V     ( I )         = V     ( J )
                  V     ( J )         = ATEMP
                  ATEMP               = BL    ( N + I )
                  BL    ( N + I )     = BL    ( N + J )
                  BL    ( N + J )     = ATEMP
                  ATEMP               = BU    ( N + I )
                  BU    ( N + I )     = BU    ( N + J )
                  BU    ( N + J )     = ATEMP
                  PNAME               = CNAME( I )
                  CNAME( I )          = CNAME( J )
                  CNAME( J )          = PNAME
                  GO TO 120
               END IF
  110       CONTINUE   
            GO TO 130
         END IF
  120 CONTINUE   
  130 CONTINUE
C
C  Add one to M for linear objective row.
C  Also set BL and BU for the objective row.
C  If the objective function has a linear part, set IOBJ to M.  
C
      M           = M + 1
      BL( N + M ) = - BIG
      BU( N + M ) = BIG
      X( N + M ) = zero
      IOBJ        = 0
      IF ( NNOBJ .LT. N ) IOBJ = M
C
C  Set up AA( I ), KA( J ) and HA( I ).  
C  AA( I ) gives the i-th element in the Jacobian.
C  KA( J ) gives starting address in AA of entries for variable J.
C  HA( I ) gives constraint index for i-th element in AA.
C
      IF ( NDEN .EQ. 1 ) THEN
C
C  Jacobian is to be stored in dense format.
C  Ensure there is sufficient room in CWK.
C
         IF ( LCWK .LT. N ) THEN
C    It used to be M + N - LCWK instead of N - LCWK.
            IF ( out .GT. 0 ) WRITE( out, 2000 )
     *           'CWK   ','LCWK  ', N - LCWK
            STOP
         END IF
         NE   = M * N
C
C  Ensure NEMAX is big enough to store the Jacobian.
C
         IF ( NEMAX .LT. NE ) THEN
            IF ( out .GT. 0 ) WRITE( out, 2000 )
     *           'AA    ','NEMAX ', NE - NEMAX
            STOP
         END IF
C
C  Use CGR to find entries in dense Jacobian.
C
         CALL CUTEST_cgr( status, N, M, X, .FALSE., M, V, 
     *              CWK, .FALSE., M, N, AA )
C
C  Set KA( J ) and HA( J ).
C
         IF ( NE .GT. 0 ) THEN
            KA( 1 ) = 1
            DO 220 J = 1, N 
               KA( J + 1 ) = KA( J ) + M
               K           = KA( J ) - 1
               DO 210 I = 1, M
                  HA( K + I ) = I
  210          CONTINUE
C
C  Copy gradient of linear part of objective function
C  into row IOBJ of Jacobian.
C
               IF ( J .GT. NNOBJ ) AA( K + IOBJ ) = CWK( J )
  220       CONTINUE
         END IF
      ELSE
C
C  Jacobian is to be stored in sparse format.
C  Ensure there is sufficient room in CWK.
C
         IF ( LCWK .LT. NEMAX ) THEN
            IF ( out .GT. 0 ) WRITE( out, 2000 )
     *           'CWK   ','LCWK  ', NEMAX - LCWK
            STOP
         END IF
C
C  Partition the integer sparse work vector ICWK.
C
         JSTRT = 0
         INDV  = JSTRT + N + 1
         INDF  = INDV  + NEMAX
C
C  Ensure there is sufficient room in ICWK.
C
         IF ( LICWK .LT. INDF + NEMAX ) THEN
            IF ( out .GT. 0 ) WRITE( out, 2000 )
     *         'ICWK  ','LICWK', INDF + NEMAX - LICWK
            STOP
         END IF
C
C  Use CSGR to find entries in sparse Jacobian.
C  Since CSGR and MINOS use different sparse formats,
C  store Jacobian temporarily in CWK.
C
         CALL CUTEST_csgr( status, N, M, .FALSE., M, V, X, NE, NEMAX,
     *              CWK, ICWK( INDV + 1 ), ICWK( INDF + 1 ) )
         K = NE
C
C  Initialize KA
C
         DO 250 J = 1, N
            KA( J ) = 0
  250    CONTINUE
C
C  Count Jacobian entries for each variable J.
C  Store counts in KA( J ).
C  Don't include nonlinear objective function entries.
C

         DO 300 II = 1, NE
            J = ICWK( INDV + II )
            I = ICWK( INDF + II )
            IF ( I .GT. 0 .OR. J .GT. NNOBJ ) THEN
               KA( J ) = KA( J ) + 1
            ELSE
               K = K - 1
            END IF
  300    CONTINUE
         KA( N + 1 ) = K + 1 
C
C  Now set KA( J ) to starting address for variable J.
C
         DO 310 J = N, 1, -1
            KA  (         J ) = KA( J + 1 ) - KA( J )
            ICWK( JSTRT + J ) = 0
  310    CONTINUE
C
C  Loop through nonlinear Jacobian entries.
C  Put correct entries in AA and HA.
C  Use KA to keep track of position for each variable J.
C  Also count nonlinear Jacobian entries for each variable J.
C  Store count in ICWK( JSTRT + J ).
C  
         NJAC = 0
         DO 320 K = 1, NE
            J = ICWK( INDV + K )
            I = ICWK( INDF + K )
            IF ( I .GT. 0 .AND. I .LE. NNCON .AND. J .LE. NNJAC ) THEN
               II                = KA ( J )
               AA  (        II ) = CWK( K )
               HA  (        II ) = I
               KA  (         J ) = II   + 1
               ICWK( JSTRT + J ) = ICWK( JSTRT + J ) + 1
               NJAC              = NJAC + 1
            END IF
  320    CONTINUE
C
C  Now loop through linear Jacobian entries,
C  including linear objective function entries.
C  Put correct entries in AA and HA.
C  Use KA to keep track of position for each variable J.
C  
         DO 330 K = 1, NE
            J = ICWK( INDV + K )
            I = ICWK( INDF + K )
            IF ( I .EQ. 0 .AND. J .GT. NNOBJ ) THEN
               II       = KA ( J )
               AA( II ) = CWK( K )
               HA( II ) = IOBJ
               KA(  J ) = II + 1
            ELSE IF ( I .GT. NNCON .OR. 
     *         ( I .GT. 0 .AND. J .GT. NNJAC ) ) THEN
               II       = KA ( J )
               AA( II ) = CWK( K )
               HA( II ) = I
               KA(  J ) = II + 1
            END IF
  330    CONTINUE
C
C  Reset KA( J ) and set ICWK( JSTRT + J ).
C  ICWK( JSTRT + J ) now gives starting address 
C  for variable J in nonlinear Jacobian.
C  These addresses are needed in FUNCON.
C
         ICWK( JSTRT + N + 1 ) = NJAC + 1
         DO 340 J = N, 2, -1
            KA  (         J ) = KA  ( J - 1 )
            ICWK( JSTRT + J ) = ICWK( JSTRT + J + 1 )
     *                          - ICWK( JSTRT + J )
  340    CONTINUE
         KA  (         1 ) = 1 
         ICWK( JSTRT + 1 ) = 1
         NE = KA ( N + 1 ) - 1
      END IF
      OBJADD = zero
      DO 400 IG = 1, NG 
         I = IWK( KNDOFC + IG )
C
C  Incorporate nonzero constants from linear constraints 
C  as bounds on slack variables.  (Constants for nonlinear
C  constraints are added in CCFG or CCFSG, which are called
C  by FUNCON.)
C
         IF ( I .GT. 0 ) THEN
            JSLACK = N + I
            IF ( I .GT. NNCON ) THEN
               IF ( WK( B + IG ) .NE. zero ) THEN
                  BU( JSLACK ) = BU( JSLACK ) - WK( B + IG ) 
     *                           * WK( GSCALE + IG )
                  BL( JSLACK ) = BL( JSLACK ) - WK( B + IG )
     *                           * WK( GSCALE + IG )
               END IF
            END IF
C
C  If possible, set slack variables to be nonbasic at zero.
C
            X ( JSLACK ) = MAX( zero, BL( JSLACK ) )
            X ( JSLACK ) = MIN( X( JSLACK ), BU( JSLACK ) )
C
C  Incorporate nonzero constants from objective function groups
C  only if objective function is completely linear.  (If objective
C  function has nonlinear part, constants are added in COFG, which
C  is called by FUNOBJ.)
C
         ELSE IF ( NNOBJ .EQ. 0 ) THEN
            OBJADD = OBJADD - WK( B + IG ) * WK( GSCALE + IG )
         END IF
  400 CONTINUE
C
C  Assign names to problem, constraints, and objective function.
C
      NAMES( 1 ) = PNAME( 1 : 8 )
      LTEMP   = .FALSE.
      DO 500 IG = 1, NG 
         I = IWK( KNDOFC + IG )
         IF ( I .LE. 0 .AND. .NOT. LTEMP ) THEN
C
C  Name objective function after its first group.
C
            NAMES( 2 )  = CHA( GNAMES + IG )
            LTEMP      = .TRUE.
         END IF
  500 CONTINUE
      NAMES( 3 )  = 'RHS     '
      NAMES( 4 )  = 'Ranges  '
      NAMES( 5 )  = 'Bounds  '
      RETURN
C
C  Non-executable statements.
C
 2000 FORMAT( /, ' ** SUBROUTINE MINOS_setup: array length ', A6, 
     *        ' too small.', /, ' -- Minimization abandoned.',
     *        /, ' -- Increase the parameter ', A6, ' by at least ', I8,
     *           ' and restart.'  )

C  End of MINOS_setup

      END

      SUBROUTINE FUNOBJ( MODE, n, X, f, G, nstate, nprob, Z, nwcore )
      INTEGER ::  mode, n, nstate, nprob, nwcore
      DOUBLE PRECISION :: f
      DOUBLE PRECISION :: X( n ), G( n ), Z( nwcore )

C  Local variables

      LOGICAL            GRAD
      IF ( MODE .EQ. 0 ) THEN
         GRAD = .FALSE.
      ELSE
         GRAD = .TRUE.
      END IF
      CALL CUTEST_cofg( status, N, X, F, G, GRAD )
      RETURN
      END

      SUBROUTINE FUNCON( MODE, M, N, NJAC, X, F, G, NSTATE, NPROB,
     *                   Z, NWCORE )
      INTEGER            MODE, M, N, NJAC, NSTATE, NPROB, NWCORE
CS    REAL               X( N ), F( M ), G( NJAC ), Z( NWCORE )
CD    DOUBLE PRECISION   X( N ), F( M ), G( NJAC ), Z( NWCORE )
C
C MINOS common block
C
      INTEGER            NCOM, NDEN, NLAG, NMAJOR, NMINOR
CS    REAL               PENPAR, ROWTOL
CD    DOUBLE PRECISION   PENPAR, ROWTOL
      COMMON    /M8AL1 / PENPAR, ROWTOL, NCOM, NDEN, NLAG,
     *                   NMAJOR, NMINOR
C
C  Sparse Jacobian common block 
C
      INTEGER            LICWK, LCWK
C    to accommodate increased MMAX for large installations.
CBIG  PARAMETER        ( LICWK = 165001 )
CMED  PARAMETER        ( LICWK = 16501 )
CTOY  PARAMETER        ( LICWK = 3301 )
CCUS  PARAMETER        ( LICWK = 165001 )
      INTEGER            ICWK( LICWK )
C    to accommodate increased MMAX for large installations.
CBIG  PARAMETER        ( LCWK  = 75000 )
CMED  PARAMETER        ( LCWK  = 7500 )
CTOY  PARAMETER        ( LCWK  = 1500 )
CCUS  PARAMETER        ( LCWK  = 75000 )
CS    REAL               CWK ( LCWK )
CD    DOUBLE PRECISION   CWK ( LCWK )
      INTEGER            JSTRT, INDV, INDF 
      COMMON / SPJAC /   CWK, ICWK, JSTRT, INDV, INDF
      SAVE             / SPJAC /
C
C Local variables
C
      INTEGER            I, J, K, NNZJ
      LOGICAL            GRAD
C
      IF ( MODE .EQ. 0 ) THEN
         GRAD = .FALSE.
      ELSE
         GRAD = .TRUE.
      END IF
C
      IF ( NDEN .EQ. 1 ) THEN
C
C  Jacobian is stored in dense format.
C
         CALL CUTEST_ccfg( status, N, M, X, M, F, .FALSE., M, N, G, GRAD )
      ELSE
C
C  Jacobian is stored in sparse format.
C
         CALL CUTEST_ccfsg( status, N, M, X, M, F, NNZJ, NJAC, CWK, 
     *               ICWK( INDV + 1 ), ICWK( INDF + 1 ), GRAD )
         IF ( GRAD ) THEN
C
C  Copy Jacobian from CCFSG, contained in CWK,
C  into MINOS Jacobian G in correct order.
C  Use ICWK( JSTRT + J ) to keep track of position for variable J.  
C
            DO 130 I = 1, NNZJ
               J                 = ICWK( INDV  + I )
               K                 = ICWK( JSTRT + J )
               G   (         K ) = CWK ( I )
               ICWK( JSTRT + J ) = K + 1
  130       CONTINUE
C
C  Reset ICWK( JSTRT + J ).
C
            DO 140 J = N, 2, -1
               ICWK( JSTRT + J ) = ICWK( JSTRT + J - 1 )
  140       CONTINUE
            ICWK( JSTRT + 1 ) = 1
         END IF
      END IF
C
      RETURN
      END

      SUBROUTINE MATMOD( NCYCLE, NPROB, FINISH, M, N, NB, NE, NKA, NS, 
     *                   NSCL, A, HA, KA, BL, BU, ASCALE, HS, ID1, ID2,
     *                   X, PI, Z, NWCORE )
      INTEGER :: ncycle, nprob, m, n, nb, ne, nka, ns, nscl, nwcore
      INTEGER ::  HA( ne ), HS( nb )
      INTEGER ::  KA( nka ), ID1( nb ), ID2( nb )
      DOUBLE PRECISION :: A( ne ), ASCALE( nscl ), BL( nb ), BU( nb ),
     *                   X( nb ), PI( m ), Z( nwcore )
      LOGICAL :: finish

C  MINOS common block

      INTEGER :: iread, iprint, isumm
      COMMON    /M1FILE/ iread, iprint, isumm
      IF ( iprint .GT. 0 ) WRITE( iprint, 2000 )
      IF ( isumm .GT. 0 ) WRITE( isumm, 2000 )
      finish = .TRUE.
      RETURN
C
C  Non-executable statements
C
 2000 FORMAT(/ 'Subroutine MATMOD has not been loaded.')
      END
