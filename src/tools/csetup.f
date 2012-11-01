! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CSETUP( INPUT, IOUT, N, M, X, BL, BU, NMAX, &
                         EQUATN, LINEAR, V, CL, CU, MMAX, EFIRST, &
                         LFIRST, NVFRST )
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: INPUT, IOUT, N, M, NMAX, MMAX
      LOGICAL :: EFIRST, LFIRST, NVFRST
      REAL ( KIND = wp ) :: X ( NMAX ), BL ( NMAX ), BU ( NMAX )
      REAL ( KIND = wp ) :: V ( MMAX ), CL ( MMAX ), CU ( MMAX )
      LOGICAL :: EQUATN( MMAX ), LINEAR( MMAX )

!  Set up the input data for the constrained optimization tools.

!  Nick Gould, for CGT productions,
!  30th October, 1991.
!  Ingrid Bongartz added option to reorder variables, so that
!  nonlinear variables come first.  Within the nonlinear variables,
!  the smaller set of either the nonlinear objective or nonlinear
!  Jacobian variables appears first.
!  5th November, 1992.

!  Workspace arrays.

      INTEGER :: LIWK, LWK, LFUVAL, LLOGIC, LCHARA

! ---------------------------------------------------------------------

!  Parameters whose value might be changed by the user:

!  The following parameters define the sizes of problem
!  dependent arrays. These may be changed by the user to
!  suit a particular problem or system configuration.

!  The TOOLS will issue error messages if any of these sizes
!  is too small, telling which parameter to increase.

! ---------------------------------------------------------------------

      INCLUDE 'tools.siz'

      INTEGER :: IWK( LIWK )
      LOGICAL :: LOGI ( LLOGIC )
      CHARACTER ( LEN = 10 ) :: CHA ( LCHARA )
      REAL ( KIND = wp ) :: WK ( LWK )
      REAL ( KIND = wp ) :: FUVALS ( LFUVAL )

! ---------------------------------------------------------------------

!  End of parameters which might be changed by the user.

! ---------------------------------------------------------------------

!  Integer variables from the GLOBAL common block.

      INTEGER :: NG, NELNUM, NGEL, NVARS, NNZA, NGPVLU
      INTEGER :: NEPVLU, NG1, NEL1, ISTADG, ISTGP, ISTADA
      INTEGER :: ISTAEV, ISTEP, ITYPEG, KNDOFC, ITYPEE
      INTEGER :: IELING, IELVAR, ICNA, ISTADH, INTVAR, IVAR
      INTEGER :: ICALCF, ITYPEV, IWRK, A, B
      INTEGER :: U, GPVALU, EPVALU
      INTEGER :: ESCALE, GSCALE, VSCALE, GVALS, XT, DGRAD
      INTEGER :: Q, WRK, INTREP, GXEQX, GNAMES, VNAMES
      INTEGER :: LO, CH, LIWORK, LWORK, NGNG, FT
      INTEGER :: LA, LB, NOBJGR, LU, LELVAR
      INTEGER :: LSTAEV, LSTADH, LNTVAR, LCALCF
      INTEGER :: LELING, LINTRE, LFT, LGXEQX, LSTADG, LGVALS
      INTEGER :: LICNA, LSTADA, LKNDOF, LGPVLU, LEPVLU
      INTEGER :: LGSCAL, LESCAL, LVSCAL, LCALCG

!  Integer variables from the LOCAL common block.

      INTEGER :: LFXI, LGXI, LHXI, LGGFX, LDX, LGRJAC
      INTEGER :: LQGRAD, LBREAK, LP, LXCP, LX0, LGX0
      INTEGER :: LDELTX, LBND, LWKSTR, LSPTRS, LSELTS, LINDEX
      INTEGER :: LSWKSP, LSTAGV, LSTAJC, LIUSED, LFREEC
      INTEGER :: LNNONZ, LNONZ2, LSYMMD, LSYMMH
      INTEGER :: LSLGRP, LSVGRP, LGCOLJ, LVALJR, LSEND
      INTEGER :: LNPTRS, LNELTS, LNNDEX, LNWKSP, LNSTGV
      INTEGER :: LNSTJC, LNIUSE, LNFREC, LNNNON, LNNNO2, LNSYMD
      INTEGER :: LNSYMH, LNLGRP, LNVGRP, LNGCLJ, LNVLJR, LNQGRD
      INTEGER :: LNBRAK, LNP, LNBND, LNFXI, LNGXI, LNGUVL
      INTEGER :: LNHXI, LNHUVL, LNGGFX, LNDX, LNGRJC, LIWK2
      INTEGER :: LWK2, MAXSIN, NINVAR, MAXSEL, LNIWTR
      INTEGER :: NTYPE, NSETS, LSTYPE, LSSWTR, LSSIWT
      INTEGER :: LSIWTR, LSWTRA, LNTYPE, LNSWTR, LNSIWT
      INTEGER :: LNWTRA, LSISET, LSSVSE, LNISET, LNSVSE
      LOGICAL :: ALTRIV, FIRSTG

!  variables from the PRFCTS common block

      INTEGER :: NC2OF, NC2OG, NC2OH,  NC2CF,  NC2CG,  NC2CH
      INTEGER :: NHVPR, PNC
      REAL :: SUTIME, STTIME
      COMMON / GLOBAL /  IWK, WK, FUVALS, LOGI, &
                         NG, NELNUM, NGEL, NVARS, NNZA, NGPVLU, &
                         NEPVLU, NG1, NEL1, ISTADG, ISTGP, ISTADA, &
                         ISTAEV, ISTEP, ITYPEG, KNDOFC, ITYPEE, &
                         IELING, IELVAR, ICNA, ISTADH, INTVAR, IVAR, &
                         ICALCF, ITYPEV, IWRK, A, B, &
                         U, GPVALU, EPVALU, &
                         ESCALE, GSCALE, VSCALE, GVALS, XT, DGRAD, &
                         Q, WRK, INTREP, GXEQX, GNAMES, VNAMES, &
                         LO, CH, LIWORK, LWORK, NGNG, FT, &
                         ALTRIV, FIRSTG, &
                         LA, LB, NOBJGR, LU, LELVAR, &
                         LSTAEV, LSTADH, LNTVAR, LCALCF, &
                         LELING, LINTRE, LFT, LGXEQX, LSTADG, LGVALS, &
                         LICNA, LSTADA, LKNDOF, LGPVLU, LEPVLU, &
                         LGSCAL, LESCAL, LVSCAL, LCALCG
      COMMON / CHARA /   CHA
      COMMON / LOCAL /   LFXI, LGXI, LHXI, LGGFX, LDX, LGRJAC, &
                         LQGRAD, LBREAK, LP, LXCP, LX0, LGX0, &
                         LDELTX, LBND, LWKSTR, LSPTRS, LSELTS, LINDEX, &
                         LSWKSP, LSTAGV, LSTAJC, LIUSED, LFREEC, &
                         LNNONZ, LNONZ2, LSYMMD, LSYMMH, &
                         LSLGRP, LSVGRP, LGCOLJ, LVALJR, LSEND, &
                         LNPTRS, LNELTS, LNNDEX, LNWKSP, LNSTGV, &
                         LNSTJC, LNIUSE, LNFREC, LNNNON, LNNNO2, LNSYMD, &
                         LNSYMH, LNLGRP, LNVGRP, LNGCLJ, LNVLJR, LNQGRD, &
                         LNBRAK, LNP, LNBND, LNFXI, LNGXI, LNGUVL, &
                         LNHXI, LNHUVL, LNGGFX, LNDX, LNGRJC, LIWK2, &
                         LWK2, MAXSIN, NINVAR, MAXSEL, NTYPE, &
                         NSETS, LSTYPE, LSSWTR, LSSIWT, LSIWTR, &
                         LSWTRA, LNTYPE, LNSWTR, LNSIWT, LNIWTR, &
                         LNWTRA, LSISET, LSSVSE, LNISET, LNSVSE
      INTEGER :: NNOV, NNJV 
      COMMON / NNVARS /  NNOV, NNJV 
      INTEGER :: IOUT2
      COMMON / OUTPUT /  IOUT2
      COMMON / PRFCTS /  NC2OF, NC2OG, NC2OH, NC2CF, NC2CG, NC2CH, &
                         NHVPR, PNC, SUTIME, STTIME
      INTEGER :: NUMVAR, NUMCON
      COMMON / DIMS /    NUMVAR, NUMCON
      SAVE             / GLOBAL /, / LOCAL /, / CHARA /, / NNVARS /, &
                       / OUTPUT /, / DIMS /, / PRFCTS /

!  Local variables.

      INTEGER :: IALGOR, IPRINT, INFORM, I, IG, J, JG, MEND
      INTEGER :: MEQ, MLIN, NSLACK, NELTYP, NGRTYP
      INTEGER :: II, K, IEL, JWRK, KNDV, NNLIN, NEND
      INTEGER :: ITEMP
      LOGICAL :: FDGRAD, DEBUG, LTEMP
      REAL :: DUM,    CPUTIM
      CHARACTER ( LEN = 8 ) :: PNAME
      CHARACTER ( LEN = 10 ) :: CTEMP
      CHARACTER ( LEN = 10 ) :: CHTEMP
      REAL ( KIND = wp ) :: ATEMP, ZERO
      PARAMETER ( ZERO = 0.0_wp )
      REAL ( KIND = wp ) :: OBFBND( 2 )
      EXTERNAL :: RANGE, CPUTIM
      SUTIME = CPUTIM( DUM )
      IOUT2 = IOUT
      DEBUG = .FALSE.
      DEBUG = DEBUG .AND. IOUT > 0
      IPRINT = 0
      IF ( DEBUG ) IPRINT = 3

!  Input the problem dimensions.

      READ( INPUT, 1001 ) N, NG, NELNUM, NGEL, NVARS, NNZA, NGPVLU, &
                          NEPVLU, NELTYP, NGRTYP
      IF ( N <= 0 ) THEN
         CLOSE( INPUT )
         IF ( IOUT > 0 ) WRITE( IOUT, 2020 )
         STOP
      END IF
      IF ( NG <= 0 ) THEN
         CLOSE( INPUT )
         IF ( IOUT > 0 ) WRITE( IOUT, 2030 )
         STOP
      END IF
      IF ( N > NMAX ) THEN
         CLOSE( INPUT )
         IF ( IOUT > 0 ) THEN
            WRITE( IOUT, 2000 ) 'X   ', 'NMAX  ', N - NMAX
            WRITE( IOUT, 2000 ) 'BL  ', 'NMAX  ', N - NMAX
            WRITE( IOUT, 2000 ) 'BU  ', 'NMAX  ', N - NMAX
         END IF
         STOP
      END IF

!  Input the problem type.

      READ( INPUT, 1000 ) IALGOR, PNAME

!  Set useful integer values.

      NG1 = NG + 1
      NGNG = NG + NG
      NEL1 = NELNUM + 1

!  Partition the integer workspace.

      ISTADG = 0
      ISTGP = ISTADG + NG1
      ISTADA = ISTGP + NG1
      ISTAEV = ISTADA + NG1
      ISTEP = ISTAEV + NEL1
      ITYPEG = ISTEP + NEL1
      KNDOFC = ITYPEG + NG
      ITYPEE = KNDOFC + NG
      IELING = ITYPEE + NELNUM
      IELVAR = IELING + NGEL
      ICNA = IELVAR + NVARS
      ISTADH = ICNA + NNZA
      INTVAR = ISTADH + NEL1
      IVAR = INTVAR + NEL1
      ICALCF = IVAR + N
      ITYPEV = ICALCF + MAX( NELNUM, NG )
      IWRK = ITYPEV + N
      LIWORK = LIWK - IWRK

!  Ensure there is sufficient room.

      IF ( LIWORK < 0 ) THEN
         CLOSE( INPUT )
         IF ( IOUT > 0 ) WRITE( IOUT, 2000 ) &
             'IWK   ', 'LIWK  ', - LIWORK
         STOP
      END IF

!  Partition the real workspace.

      A = 0
      B = A + NNZA
      U = B + NG
      GPVALU = U + NG
      EPVALU = GPVALU + NGPVLU
      ESCALE = EPVALU + NEPVLU
      GSCALE = ESCALE + NGEL
      VSCALE = GSCALE + NG
      GVALS = VSCALE + N
      XT = GVALS + 3 * NG
      DGRAD = XT + N
      Q = DGRAD + N
      FT = Q + N
      WRK = FT + NG
      LWORK = LWK - WRK

!  Ensure there is sufficient room.

      IF ( LWORK < 0 ) THEN
         CLOSE( INPUT )
         IF ( IOUT > 0 ) WRITE( IOUT, 2000 ) &
             'WK   ', 'LWK   ', - LWORK
         STOP
      END IF

!  Partition the logical workspace.

      INTREP = 0
      GXEQX = INTREP + NELNUM
      LO = GXEQX + NGNG

!  Ensure there is sufficient room.

      IF ( LLOGIC < LO ) THEN
         CLOSE( INPUT )
         IF ( IOUT > 0 ) WRITE( IOUT, 2000 ) &
             'LOGI  ', 'LLOGIC', LO - LLOGIC
         STOP
      END IF

!  Partition the character workspace.

      GNAMES = 0
      VNAMES = GNAMES + NG
      CH = VNAMES + N

!  Ensure there is sufficient room.

      IF ( LCHARA < CH + 1 ) THEN
         CLOSE( INPUT )
         IF ( IOUT > 0 ) WRITE( IOUT, 2000 ) &
             'CHA   ', 'LCHARA', CH + 1 - LCHARA
         STOP
      END IF

!  Record the lengths of arrays.

      LSTADG = MAX( 1, NG1 )
      LSTADA = MAX( 1, NG1 )
      LSTAEV = MAX( 1, NEL1 )
      LKNDOF = MAX( 1, NG )
      LELING = MAX( 1, NGEL )
      LELVAR = MAX( 1, NVARS )
      LICNA = MAX( 1, NNZA )
      LSTADH = MAX( 1, NEL1 )
      LNTVAR = MAX( 1, NEL1 )
      LCALCF = MAX( 1, NELNUM, NG )
      LCALCG = MAX( 1, NG )
      LA = MAX( 1, NNZA )
      LB = MAX( 1, NG )
      LU = MAX( 1, NG )
      LESCAL = MAX( 1, NGEL )
      LGSCAL = MAX( 1, NG )
      LVSCAL = MAX( 1, N )
      LFT = MAX( 1, NG )
      LGVALS = MAX( 1, NG )
      LINTRE = MAX( 1, NELNUM )
      LGXEQX = MAX( 1, NGNG )
      LGPVLU = MAX( 1, NGPVLU )
      LEPVLU = MAX( 1, NEPVLU )
!     LSTGP = MAX( 1, NG1 )
!     LSTEP = MAX( 1, NEL1 )
!     LTYPEG = MAX( 1, NG )
!     LTYPEE = MAX( 1, NELNUM )
!     LIVAR = MAX( 1, N )
!     LBL = MAX( 1, N )
!     LBU = MAX( 1, N )
!     LX = MAX( 1, N )
!     LXT = MAX( 1, N )
!     LDGRAD = MAX( 1, N )
!     LQ = MAX( 1, N )


!  Print out problem data. input the number of variables, groups,
!  elements and the identity of the objective function group.

      IF ( IALGOR == 2 ) THEN
         READ( INPUT, 1002 ) NSLACK, NOBJGR
      ELSE
         NSLACK = 0
      END IF
      IF ( DEBUG ) WRITE( IOUT, 1100 ) PNAME, N, NG, NELNUM
      CHA( CH + 1 ) = PNAME // '  '

!  Input the starting addresses of the elements in each group,
!  of the parameters used for each group and
!  of the nonzeros of the linear element in each group.

      READ( INPUT, 1010 ) ( IWK( ISTADG + I ), I = 1, NG1 )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'ISTADG', &
 ( IWK( ISTADG + I ), I = 1, NG1 )
      READ( INPUT, 1010 ) ( IWK( ISTGP + I ), I = 1, NG1 )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'ISTGP ', &
 ( IWK( ISTGP + I ), I = 1, NG1 )
      READ( INPUT, 1010 ) ( IWK( ISTADA + I ), I = 1, NG1 )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'ISTADA', &
 ( IWK( ISTADA + I ), I = 1, NG1 )

!  Input the starting addresses of the variables and parameters
!  in each element.

      READ( INPUT, 1010 ) ( IWK( ISTAEV + I ), I = 1, NEL1 )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'ISTAEV', &
 ( IWK( ISTAEV + I ), I = 1, NEL1 )
      READ( INPUT, 1010 ) ( IWK( ISTEP + I ), I = 1, NEL1 )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'ISTEP ', &
 ( IWK( ISTEP + I ), I = 1, NEL1 )

!  Input the group type of each group

      READ( INPUT, 1010 ) ( IWK( ITYPEG + I ), I = 1, NG )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'ITYPEG', &
 ( IWK( ITYPEG + I ), I = 1, NG )
      IF ( IALGOR >= 2 ) THEN
         READ( INPUT, 1010 ) ( IWK( KNDOFC + I ), I = 1, NG )
         IF ( DEBUG ) WRITE( IOUT, 1110 ) 'KNDOFC', &
 ( IWK( KNDOFC + I ), I = 1, NG )
      END IF

!  Input the element type of each element

      READ( INPUT, 1010 ) ( IWK( ITYPEE + I ), I = 1, NELNUM )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'ITYPEE', &
 ( IWK( ITYPEE + I ), I = 1, NELNUM )

!  Input the number of internal variables for each element.

      READ( INPUT, 1010 ) ( IWK( INTVAR + I ), I = 1, NELNUM )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'INTVAR', &
 ( IWK( INTVAR + I ), I = 1, NELNUM )

!  Input the identity of each individual element.

      READ( INPUT, 1010 ) ( IWK( IELING + I ), I = 1, NGEL )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'IELING', &
 ( IWK( IELING + I ), I = 1, NGEL )

!  Input the variables in each group's elements.

      NVARS = IWK( ISTAEV + NEL1 ) - 1
      READ( INPUT, 1010 ) ( IWK( IELVAR + I ), I = 1, NVARS )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'IELVAR', &
 ( IWK( IELVAR + I ), I = 1, NVARS )

!  Input the column addresses of the nonzeros in each linear element.

      READ( INPUT, 1010 ) ( IWK( ICNA + I ), I = 1, NNZA )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'ICNA  ', &
 ( IWK( ICNA + I ), I = 1, NNZA )

!  Input the values of the nonzeros in each linear element, the
!  constant term in each group, the lower and upper bounds on
!  the variables.

      READ( INPUT, 1020 ) ( WK( A + I ), I = 1, NNZA )
      IF ( DEBUG ) WRITE( IOUT, 1120 ) 'A     ', &
 ( WK( A + I ), I = 1, NNZA )
      READ( INPUT, 1020 ) ( WK( B + I ), I = 1, NG )
      IF ( DEBUG ) WRITE( IOUT, 1120 ) 'B     ', &
 ( WK( B + I ), I = 1, NG )
      IF ( IALGOR <= 2 ) THEN
         READ( INPUT, 1020 ) ( BL( I ), I = 1, N )
         IF ( DEBUG ) WRITE( IOUT, 1120 ) 'BL    ', &
 ( BL( I ), I = 1, N )
         READ( INPUT, 1020 ) ( BU( I ), I = 1, N )
         IF ( DEBUG ) WRITE( IOUT, 1120 ) 'BU    ', &
 ( BU( I ), I = 1, N )
      ELSE

!  Use GVALS and FT as temporary storage for the constraint bounds.

         READ( INPUT, 1020 ) ( BL( I ), I = 1, N ), &
 ( WK( GVALS + I ), I = 1, NG )
         IF ( DEBUG ) WRITE( IOUT, 1120 ) 'BL    ', &
 ( BL( I ), I = 1, N ), ( WK( GVALS + I ), I = 1, NG )
         READ( INPUT, 1020 ) ( BU( I ), I = 1, N ), &
 ( WK( FT + I ), I = 1, NG )
         IF ( DEBUG ) WRITE( IOUT, 1120 ) 'BU    ', &
 ( BU( I ), I = 1, N ), ( WK( FT + I ), I = 1, NG )
      END IF

!   Input the starting point for the minimization.

      READ( INPUT, 1020 ) ( X( I ), I = 1, N )
      IF ( DEBUG ) WRITE( IOUT, 1120 ) 'X     ', &
 ( X( I ), I = 1, N )
      IF ( IALGOR >= 2 ) THEN
         READ( INPUT, 1020 )( WK( U + I ), I = 1, NG )
         IF ( DEBUG ) WRITE( IOUT, 1120 ) 'U     ', &
 ( WK( U + I ), I = 1, NG )
      END IF

!  Input the parameters in each group.

      READ( INPUT, 1020 ) ( WK( GPVALU + I ), I = 1, NGPVLU )
      IF ( DEBUG ) WRITE( IOUT, 1120 ) 'GPVALU', &
 ( WK( GPVALU + I ), I = 1, NGPVLU )

!  Input the parameters in each individual element.

      READ( INPUT, 1020 ) ( WK( EPVALU + I ), I = 1, NEPVLU )
      IF ( DEBUG ) WRITE( IOUT, 1120 ) 'EPVALU', &
 ( WK( EPVALU + I ), I = 1, NEPVLU )

!  Input the scale factors for the nonlinear elements.

      READ( INPUT, 1020 ) ( WK( ESCALE + I ), I = 1, NGEL )
      IF ( DEBUG ) WRITE( IOUT, 1120 ) 'ESCALE', &
 ( WK( ESCALE + I ), I = 1, NGEL )

!  Input the scale factors for the groups.

      READ( INPUT, 1020 ) ( WK( GSCALE + I ), I = 1, NG )
      IF ( DEBUG ) WRITE( IOUT, 1120 ) 'GSCALE', &
 ( WK( GSCALE + I ), I = 1, NG )

!  Input the scale factors for the variables.

      READ( INPUT, 1020 ) ( WK( VSCALE + I ), I = 1, N )
      IF ( DEBUG ) WRITE( IOUT, 1120 ) 'VSCALE', &
 ( WK( VSCALE + I ), I = 1, N )

!  Input the lower and upper bounds on the objective function.

      READ( INPUT, 1080 ) OBFBND( 1 ), OBFBND( 2 )
      IF ( DEBUG ) WRITE( IOUT, 1180 ) 'OBFBND', &
          OBFBND( 1 ), OBFBND( 2 )

!  Input a logical array which says whether an element has internal
!  variables.

      READ( INPUT, 1030 ) ( LOGI( INTREP + I ), I = 1, NELNUM )
      IF ( DEBUG ) WRITE( IOUT, 1130 ) 'INTREP', &
 ( LOGI( INTREP + I ), I = 1, NELNUM )

!  Input a logical array which says whether a group is trivial.

      READ( INPUT, 1030 ) ( LOGI( GXEQX + I ), I = 1, NG )
      IF ( DEBUG ) WRITE( IOUT, 1130 ) 'GXEQX ', &
 ( LOGI( GXEQX + I ), I = 1, NG )

!  Input the names given to the groups and to the variables.

      READ( INPUT, 1040 ) ( CHA( GNAMES + I ), I = 1, NG )
      IF ( DEBUG ) WRITE( IOUT, 1140 ) 'GNAMES', &
 ( CHA( GNAMES + I ), I = 1, NG )
      READ( INPUT, 1040 ) ( CHA( VNAMES + I ), I = 1, N )
      IF ( DEBUG ) WRITE( IOUT, 1140 ) 'VNAMES', &
 ( CHA( VNAMES + I ), I = 1, N )

!  Dummy input for the names given to the element and group types.

      READ( INPUT, 1040 ) ( CHTEMP, I = 1, NELTYP )
      READ( INPUT, 1040 ) ( CHTEMP, I = 1, NGRTYP )

!  Input the type of each variable.

      READ( INPUT, 1010 ) ( IWK( ITYPEV + I ), I = 1, N )

!  Consider which groups are constraints. Of these, decide which are
!  equations, which are linear, allocate starting values for the
!  Lagrange multipliers and set lower and upper bounds on any
!  inequality constraints. Reset KNDOFC to point to the list of
!  constraint groups.

      M = 0
      DO 10 I = 1, NG
         IF ( IWK( KNDOFC + I ) == 1 ) THEN
            IWK( KNDOFC + I ) = 0
         ELSE
            M = M + 1
            IF ( M <= MMAX ) THEN
               V ( M ) = WK( U + I )
               LINEAR( M ) = LOGI( GXEQX + I ) .AND. IWK( ISTADG + I ) &
                             >= IWK( ISTADG + I + 1 )
               IF ( IWK( KNDOFC + I ) == 2 ) THEN
                  EQUATN( M ) = .TRUE.
                  CL ( M ) = ZERO
                  CU ( M ) = ZERO
               ELSE
                  EQUATN( M ) = .FALSE.
                  CL ( M ) = WK( GVALS + I )
                  CU ( M ) = WK( FT + I )
               END IF
            END IF
            IWK( KNDOFC + I ) = M
         END IF
   10 CONTINUE
      IF ( M == 0 .AND. IOUT > 0 ) WRITE( IOUT, 2010 )
      IF ( M > MMAX ) THEN
         CLOSE( INPUT )
         IF ( IOUT > 0 ) THEN
            WRITE( IOUT, 2000 ) 'V     ', 'MMAX  ', M - MMAX
            WRITE( IOUT, 2000 ) 'CL    ', 'MMAX  ', M - MMAX
            WRITE( IOUT, 2000 ) 'CU    ', 'MMAX  ', M - MMAX
            WRITE( IOUT, 2000 ) 'EQUATN', 'MMAX  ', M - MMAX
            WRITE( IOUT, 2000 ) 'LINEAR', 'MMAX  ', M - MMAX
         END IF
         STOP
      END IF
!                            the number of variables and constraints, resp.
      NUMVAR = N
      NUMCON = M
      IF ( NVFRST ) THEN

!  Ensure there is sufficient room in IWK to reorder variables.

         IF ( LIWORK < 2*N ) THEN
            CLOSE( INPUT )
            IF ( IOUT > 0 ) WRITE( IOUT, 2000 ) &
                'IWK   ', 'LIWK  ', LIWORK - 2*N
            STOP
         END IF
         KNDV = IWRK + 1 
         JWRK = KNDV + N

!  Initialize JWRK and KNDV.

         DO 20 J = 1, N
            IWK( KNDV + J ) = 0
            IWK( JWRK + J ) = J
   20    CONTINUE

!  Now identify and count nonlinear variables.
!  Keep separate counts for nonlinear objective and Jacobian variables.
!  IWK( KNDV + J ) = 0 ==> J linear everywhere
!  IWK( KNDV + J ) = 1 ==> J linear in objective, nonlinear in constraints
!  IWK( KNDV + J ) = 2 ==> J linear in constraints, nonlinear in objective
!  IWK( KNDV + J ) = 3 ==> J nonlinear everywhere

         NNLIN = 0
         NNOV = 0
         NNJV = 0
         DO 60 IG = 1, NG
            I = IWK( KNDOFC + IG )
            DO 40 II = IWK( ISTADG + IG ), IWK( ISTADG + IG + 1 ) - 1
               IEL = IWK( IELING + II )
               DO 30 K = IWK( ISTAEV + IEL ), IWK( ISTAEV + IEL + 1) - 1
                  J = IWK( IELVAR + K )
                  IF ( I > 0 ) THEN
                     IF ( IWK( KNDV + J ) == 0 ) THEN
                        IWK( KNDV + J ) = 1
                        NNJV = NNJV + 1
                        NNLIN = NNLIN + 1
                     ELSE IF ( IWK( KNDV + J ) == 2 ) THEN
                        IWK( KNDV + J ) = 3
                        NNJV = NNJV + 1
                     END IF
                  ELSE
                     IF ( IWK( KNDV + J ) == 0 ) THEN
                        IWK( KNDV + J ) = 2
                        NNOV = NNOV + 1
                        NNLIN = NNLIN + 1
                     ELSE IF ( IWK( KNDV + J ) == 1 ) THEN
                        IWK( KNDV + J ) = 3
                        NNOV = NNOV + 1
                     END IF
                  END IF
   30          CONTINUE
   40       CONTINUE
            IF ( .NOT. LOGI( GXEQX + IG ) ) THEN
               DO 50 II = IWK( ISTADA + IG ), IWK( ISTADA + IG + 1 ) - 1
                  J = IWK( ICNA + II )
                  IF ( I > 0 ) THEN
                     IF ( IWK( KNDV + J ) == 0 ) THEN
                        IWK( KNDV + J ) = 1
                        NNJV = NNJV + 1
                        NNLIN = NNLIN + 1
                     ELSE IF ( IWK( KNDV + J ) == 2 ) THEN
                        IWK( KNDV + J ) = 3
                        NNJV = NNJV + 1
                     END IF
                  ELSE
                     IF ( IWK( KNDV + J ) == 0 ) THEN
                        IWK( KNDV + J ) = 2
                        NNOV = NNOV + 1
                        NNLIN = NNLIN + 1
                     ELSE IF ( IWK( KNDV + J ) == 1 ) THEN
                        IWK( KNDV + J ) = 3
                        NNOV = NNOV + 1
                     END IF
                  END IF
   50          CONTINUE
            END IF
   60    CONTINUE
         IF ( NNLIN == 0 .OR. ( NNOV == N .AND. NNJV == N ) ) &
            GO TO 600
         IF ( NNLIN == N ) GO TO 500

!  Reorder the variables so that all nonlinear variables occur before
!  the linear ones.

         NEND = N

!  Run forward through the variables until a linear variable
!  is encountered.

         DO 420 I = 1, N
            IF ( I > NEND ) GO TO 430
            IF ( IWK( KNDV + I ) == 0 ) THEN

!  Variable I is linear. Now, run backwards through the
!  variables until a nonlinear one is encountered.

               DO 410 J = NEND, I, - 1
                  IF ( IWK( KNDV + J ) > 0 ) THEN 
                     NEND = J - 1

!  Interchange the data for variables I and J.

                     ITEMP = IWK ( JWRK + I )
                     IWK ( JWRK + I ) = IWK ( JWRK + J )
                     IWK ( JWRK + J ) = ITEMP
                     ITEMP = IWK ( KNDV + I )
                     IWK ( KNDV + I ) = IWK ( KNDV + J )
                     IWK ( KNDV + J ) = ITEMP
                     ATEMP = BL ( I )
                     BL ( I ) = BL ( J )
                     BL ( J ) = ATEMP
                     ATEMP = BU ( I )
                     BU ( I ) = BU ( J )
                     BU ( J ) = ATEMP
                     ATEMP = X ( I )
                     X ( I ) = X ( J )
                     X ( J ) = ATEMP
                     ATEMP = WK ( VSCALE + I )
                     WK ( VSCALE + I ) = WK ( VSCALE + J )
                     WK ( VSCALE + J ) = ATEMP
                     CTEMP = CHA ( VNAMES + I )
                     CHA ( VNAMES + I ) = CHA ( VNAMES + J )
                     CHA ( VNAMES + J ) = CTEMP
                     GO TO 420
                  END IF
  410          CONTINUE
               GO TO 430
            END IF
  420    CONTINUE
  430    CONTINUE 

!  Change entries in IELVAR and ICNA to reflect reordering of variables.

         DO 440 I = 1, NVARS
            J = IWK( IELVAR + I )
            IWK( IELVAR + I ) = IWK( JWRK + J ) 
  440    CONTINUE
         DO 450 I = 1, NNZA
            J = IWK( ICNA + I )
            IWK( ICNA + I ) = IWK( JWRK + J )
  450    CONTINUE
         DO 460 J = 1, N
            IWK( JWRK + J ) = J
  460    CONTINUE
  500    CONTINUE
         IF ( ( NNOV == NNLIN .AND. NNJV == NNLIN )  &
            .OR. ( NNOV == 0 ) .OR. ( NNJV == 0 ) ) GO TO 600

!  Reorder the nonlinear variables so that the smaller set (nonlinear
!  objective or nonlinear Jacobian) occurs at the beginning of the 
!  larger set.

         NEND = NNLIN
         IF ( NNJV <= NNOV ) THEN

!  Put the nonlinear Jacobian variables first.
!  Reset NNOV to indicate all nonlinear variables are treated as
!  nonlinear objective variables.

            NNOV = NNLIN
            DO 520 I = 1, NNLIN 
               IF ( I > NEND ) GO TO 530
               IF ( IWK( KNDV + I ) == 2 ) THEN

!  Variable I is linear in the Jacobian. Now, run backwards through the 
!  variables until a nonlinear Jacobian variable is encountered.

                  DO 510 J = NEND, I, - 1
                     IF ( IWK( KNDV + J ) == 1  &
                        .OR. IWK( KNDV + J ) == 3 ) THEN 
                        NEND = J - 1

!  Interchange the data for variables I and J.

                        ITEMP = IWK ( JWRK + I )
                        IWK ( JWRK + I ) = IWK ( JWRK + J )
                        IWK ( JWRK + J ) = ITEMP
                        ITEMP = IWK ( KNDV + I )
                        IWK ( KNDV + I ) = IWK ( KNDV + J )
                        IWK ( KNDV + J ) = ITEMP
                        ATEMP = BL ( I )
                        BL ( I ) = BL ( J )
                        BL ( J ) = ATEMP
                        ATEMP = BU ( I )
                        BU ( I ) = BU ( J )
                        BU ( J ) = ATEMP
                        ATEMP = X ( I )
                        X ( I ) = X ( J )
                        X ( J ) = ATEMP
                        ATEMP = WK ( VSCALE + I )
                        WK ( VSCALE + I ) = WK ( VSCALE + J )
                        WK ( VSCALE + J ) = ATEMP
                        CTEMP = CHA ( VNAMES + I )
                        CHA ( VNAMES + I ) = CHA ( VNAMES + J )
                        CHA ( VNAMES + J ) = CTEMP
                        GO TO 520
                     END IF
  510             CONTINUE
                  GO TO 530
               END IF
  520       CONTINUE
  530       CONTINUE
         ELSE

!  Put the nonlinear objective variables first.
!  Reset NNJV to indicate all nonlinear variables are treated as
!  nonlinear Jacobian variables.

            NNJV = NNLIN
            DO 550 I = 1, NNLIN 
               IF ( I > NEND ) GO TO 560
               IF ( IWK( KNDV + I ) == 1 ) THEN

!  Variable I is linear in the objective. Now, run backwards through the 
!  variables until a nonlinear objective variable is encountered.

                  DO 540 J = NEND, I, - 1
                     IF ( IWK( KNDV + J ) > 1 ) THEN 
                        NEND = J - 1

!  Interchange the data for variables I and J.

                        ITEMP = IWK ( JWRK + I )
                        IWK ( JWRK + I ) = IWK ( JWRK + J )
                        IWK ( JWRK + J ) = ITEMP
                        ITEMP = IWK ( KNDV + I )
                        IWK ( KNDV + I ) = IWK ( KNDV + J )
                        IWK ( KNDV + J ) = ITEMP
                        ATEMP = BL ( I )
                        BL ( I ) = BL ( J )
                        BL ( J ) = ATEMP
                        ATEMP = BU ( I )
                        BU ( I ) = BU ( J )
                        BU ( J ) = ATEMP
                        ATEMP = X ( I )
                        X ( I ) = X ( J )
                        X ( J ) = ATEMP
                        ATEMP = WK ( VSCALE + I )
                        WK ( VSCALE + I ) = WK ( VSCALE + J )
                        WK ( VSCALE + J ) = ATEMP
                        CTEMP = CHA ( VNAMES + I )
                        CHA ( VNAMES + I ) = CHA ( VNAMES + J )
                        CHA ( VNAMES + J ) = CTEMP
                        GO TO 550
                     END IF
  540             CONTINUE
                  GO TO 560
               END IF
  550       CONTINUE
  560       CONTINUE
         END IF

!  Change entries in IELVAR and ICNA to reflect reordering of variables.

         DO 580 I = 1, NVARS
            J = IWK( IELVAR + I )
            IWK( IELVAR + I ) = IWK( JWRK + J ) 
  580    CONTINUE
         DO 590 I = 1, NNZA
            J = IWK( ICNA + I )
            IWK( ICNA + I ) = IWK( JWRK + J )
  590    CONTINUE
  600    CONTINUE
      END IF

!  Partition the workspace arrays FUVALS, IWK and WK. Initialize
!  certain portions of IWK.

      FIRSTG = .TRUE.
      FDGRAD = .FALSE.
      CALL DINITW( N, NG, NELNUM, IWK(IELING + 1), LELING, IWK(ISTADG + 1), &
          LSTADG, IWK(IELVAR + 1), LELVAR, IWK(ISTAEV + 1), LSTAEV, &
          IWK(INTVAR + 1), LNTVAR, IWK(ISTADH + 1), LSTADH, &
          IWK(ICNA + 1), LICNA, IWK(ISTADA + 1), LSTADA, &
          IWK(ITYPEE + 1), LINTRE, &
          LOGI(GXEQX + 1), LGXEQX, LOGI(INTREP + 1), LINTRE, &
          LFUVAL, ALTRIV, .TRUE., FDGRAD, LFXI,LGXI,LHXI,LGGFX, &
          LDX, LGRJAC, LQGRAD, LBREAK, LP,     LXCP, LX0,  &
          LGX0, LDELTX, LBND,   LWKSTR, LSPTRS, LSELTS, LINDEX,  &
          LSWKSP, LSTAGV, LSTAJC, LIUSED, LFREEC, LNNONZ, LNONZ2,  &
          LSYMMD, LSYMMH, LSLGRP, LSVGRP, LGCOLJ, LVALJR, LSEND,  &
          LNPTRS, LNELTS, LNNDEX, LNWKSP, LNSTGV, LNSTJC, LNIUSE,  &
          LNFREC, LNNNON, LNNNO2, LNSYMD, LNSYMH, LNLGRP, LNVGRP, &
          LNGCLJ, LNVLJR, LNQGRD, LNBRAK, LNP,    LNBND, &
          LNFXI,  LNGXI,  LNGUVL, LNHXI,  LNHUVL, LNGGFX, &
          LNDX, LNGRJC, LIWK2, LWK2, MAXSIN, NINVAR, &
          NTYPE, NSETS, MAXSEL, LSTYPE, LSSWTR, LSSIWT, &
          LSIWTR, LSWTRA, LNTYPE, LNSWTR, LNSIWT, LNIWTR, &
          LNWTRA, LSISET, LSSVSE, LNISET, LNSVSE, RANGE, &
          IWK(IWRK + 1),LIWORK,WK(WRK + 1),LWORK,IPRINT,IOUT,INFORM )
      IF ( INFORM /= 0 ) STOP

!  Shift the starting addresses for the real workspace relative to WRK.

      LQGRAD = LQGRAD + WRK
      LBREAK = LBREAK + WRK
      LP = LP + WRK
      LXCP = LXCP + WRK
      LX0 = LX0 + WRK
      LGX0 = LGX0 + WRK
      LDELTX = LDELTX + WRK
      LBND = LBND + WRK
      LSWTRA = LSWTRA + WRK
      LWKSTR = LWKSTR + WRK

!  Shift the starting addresses for the integer workspace relative
!  to IWRK.

      LSPTRS = LSPTRS + IWRK
      LSELTS = LSELTS + IWRK
      LINDEX = LINDEX + IWRK
      LSWKSP = LSWKSP + IWRK
      LSTAGV = LSTAGV + IWRK
      LSTAJC = LSTAJC + IWRK
      LIUSED = LIUSED + IWRK
      LFREEC = LFREEC + IWRK
      LNNONZ = LNNONZ + IWRK
      LNONZ2 = LNONZ2 + IWRK
      LSYMMD = LSYMMD + IWRK
      LSYMMH = LSYMMH + IWRK
      LSLGRP = LSLGRP + IWRK
      LSVGRP = LSVGRP + IWRK
      LGCOLJ = LGCOLJ + IWRK
      LVALJR = LVALJR + IWRK
      LSTYPE = LSTYPE + IWRK
      LSSWTR = LSSWTR + IWRK
      LSSIWT = LSSIWT + IWRK
      LSIWTR = LSIWTR + IWRK
      LSISET = LSISET + IWRK
      LSSVSE = LSSVSE + IWRK
      LSEND = LSEND + IWRK
      IF ( .NOT. ( EFIRST .OR. LFIRST ) ) GOTO 340
!                            to RETURN if there are no constraints.
      IF ( M == 0 ) GOTO 340

!  Record which group is associated with each constraint.

      IF ( M > LIWK2 ) THEN
         WRITE( IOUT, 2040 )
         STOP
      END IF
      MEQ = 0
      MLIN = 0
      DO 100 IG = 1, NG
         I = IWK( KNDOFC + IG )
         IF ( I > 0 ) THEN
            IWK( LSEND + I ) = IG
            IF ( EQUATN( I ) ) MEQ = MEQ + 1
            IF ( LINEAR( I ) ) MLIN = MLIN + 1
         END IF
  100 CONTINUE
      IF ( LFIRST ) THEN
         IF ( MLIN == 0 .OR. MLIN == M ) GO TO 130

!  Reorder the constraints so that the linear constraints occur before the
!  nonlinear ones.

         MEND = M

!  Run forward through the constraints until a nonlinear constraint
!  is encountered.

         DO 120 I = 1, M
            IF ( I > MEND ) GO TO 130
            IG = IWK( LSEND + I )
!              write(6,*) ' group ', IG, ' type ', I, ' equal? ',
!     *                     EQUATN( I )
            IF ( .NOT. LINEAR( I ) ) THEN

!  Constraint I is nonlinear. Now, run backwards through the
!  constraints until a linear one is encountered.

               DO 110 J = MEND, I, - 1
                  JG = IWK( LSEND + J )
!                 write(6,*) ' group ', JG, ' type ', J,
!     *                      ' linear? ', LINEAR( J )
                  IF ( LINEAR( J ) ) THEN
!                    write(6,*) ' swaping constraints ', I,
!     *                         ' and ', J
                     MEND = J - 1

!  Interchange the data for constraints I and J.

                     IWK( LSEND + I ) = JG
                     IWK( LSEND + J ) = IG
                     IWK( KNDOFC + IG ) = J
                     IWK( KNDOFC + JG ) = I
                     LTEMP = LINEAR( I )
                     LINEAR( I ) = LINEAR( J )
                     LINEAR( J ) = LTEMP
                     LTEMP = EQUATN( I )
                     EQUATN( I ) = EQUATN( J )
                     EQUATN( J ) = LTEMP
                     ATEMP = V ( I )
                     V ( I ) = V ( J )
                     V ( J ) = ATEMP
                     ATEMP = CL ( I )
                     CL ( I ) = CL ( J )
                     CL ( J ) = ATEMP
                     ATEMP = CU ( I )
                     CU ( I ) = CU ( J )
                     CU ( J ) = ATEMP
                     GO TO 120
                  END IF
  110          CONTINUE
               GO TO 130
            END IF
  120    CONTINUE
  130    CONTINUE
         IF ( EFIRST ) THEN
            IF ( MEQ == 0 .OR. MEQ == M ) GO TO 260

!  Reorder the linear constraints so that the equations occur before
!  the inequalities.

            MEND = MLIN
            DO 220 I = 1, MLIN
               IF ( I > MEND ) GO TO 230
               IG = IWK( LSEND + I )
!                 write(6,*) ' group ', IG, ' type ', I, ' equation? ',
!     *                        EQUATN( I )
               IF ( .NOT. EQUATN( I ) ) THEN

!  Constraint I is an inequality. Now, run backwards through the
!  constraints until an equation is encountered.

                  DO 210 J = MEND, I, - 1
                     JG = IWK( LSEND + J )
!                    write(6,*) ' group ', JG, ' type ', J,
!     *                         ' equation? ', EQUATN( J )
                     IF ( EQUATN( J ) ) THEN
!                       write(6,*) ' swaping constraints ', I,
!     *                            ' and ', J
                        MEND = J - 1

!  Interchange the data for constraints I and J.

                        IWK( LSEND + I ) = JG
                        IWK( LSEND + J ) = IG
                        IWK( KNDOFC + IG ) = J
                        IWK( KNDOFC + JG ) = I
                        LTEMP = LINEAR( I )
                        LINEAR( I ) = LINEAR( J )
                        LINEAR( J ) = LTEMP
                        LTEMP = EQUATN( I )
                        EQUATN( I ) = EQUATN( J )
                        EQUATN( J ) = LTEMP
                        ATEMP = V ( I )
                        V ( I ) = V ( J )
                        V ( J ) = ATEMP
                        ATEMP = CL ( I )
                        CL ( I ) = CL ( J )
                        CL ( J ) = ATEMP
                        ATEMP = CU ( I )
                        CU ( I ) = CU ( J )
                        CU ( J ) = ATEMP
                        GO TO 220
                     END IF
  210             CONTINUE
                  GO TO 230
               END IF
  220       CONTINUE
  230       CONTINUE

!  Reorder the nonlinear constraints so that the equations occur
!  before the inequalities.

            MEND = M
            DO 250 I = MLIN + 1, M
               IF ( I > MEND ) GO TO 260
               IG = IWK( LSEND + I )
!                 write(6,*) ' group ', IG, ' type ', I, ' equation? ',
!     *                        EQUATN( I )
               IF ( .NOT. EQUATN( I ) ) THEN

!  Constraint I is an inequality. Now, run backwards through the
!  constraints until an equation is encountered.

                  DO 240 J = MEND, I, - 1
                     JG = IWK( LSEND + J )
!                    write(6,*) ' group ', JG, ' type ', J,
!     *                         ' equation? ', EQUATN( J )
                     IF ( EQUATN( J ) ) THEN
!                       write(6,*) ' swaping constraints ', I,
!     *                            ' and ', J
                        MEND = J - 1

!  Interchange the data for constraints I and J.

                        IWK( LSEND + I ) = JG
                        IWK( LSEND + J ) = IG
                        IWK( KNDOFC + IG ) = J
                        IWK( KNDOFC + JG ) = I
                        LTEMP = LINEAR( I )
                        LINEAR( I ) = LINEAR( J )
                        LINEAR( J ) = LTEMP
                        LTEMP = EQUATN( I )
                        EQUATN( I ) = EQUATN( J )
                        EQUATN( J ) = LTEMP
                        ATEMP = V ( I )
                        V ( I ) = V ( J )
                        V ( J ) = ATEMP
                        ATEMP = CL ( I )
                        CL ( I ) = CL ( J )
                        CL ( J ) = ATEMP
                        ATEMP = CU ( I )
                        CU ( I ) = CU ( J )
                        CU ( J ) = ATEMP
                        GO TO 250
                     END IF
  240             CONTINUE
                  GO TO 260
               END IF
  250       CONTINUE
  260       CONTINUE
         END IF
      ELSE
         IF ( EFIRST ) THEN
            IF ( MEQ == 0 .OR. MEQ == M ) GO TO 330

!  Reorder the constraints so that the equations occur before the
!  inequalities.

            MEND = M
            DO 320 I = 1, M
               IF ( I > MEND ) GO TO 330
               IG = IWK( LSEND + I )
!                 write(6,*) ' group ', IG, ' type ', I, ' equation? ',
!     *                        EQUATN( I )
               IF ( .NOT. EQUATN( I ) ) THEN

!  Constraint I is an inequality. Now, run backwards through the
!  constraints until an equation is encountered.

                  DO 310 J = MEND, I, - 1
                     JG = IWK( LSEND + J )
!                    write(6,*) ' group ', JG, ' type ', J,
!     *                         ' equation? ', EQUATN( J )
                     IF ( EQUATN( J ) ) THEN
!                       write(6,*) ' swaping constraints ', I,
!     *                            ' and ', J
                        MEND = J - 1

!  Interchange the data for constraints I and J.

                        IWK( LSEND + I ) = JG
                        IWK( LSEND + J ) = IG
                        IWK( KNDOFC + IG ) = J
                        IWK( KNDOFC + JG ) = I
                        LTEMP = LINEAR( I )
                        LINEAR( I ) = LINEAR( J )
                        LINEAR( J ) = LTEMP
                        LTEMP = EQUATN( I )
                        EQUATN( I ) = EQUATN( J )
                        EQUATN( J ) = LTEMP
                        ATEMP = V ( I )
                        V ( I ) = V ( J )
                        V ( J ) = ATEMP
                        ATEMP = CL ( I )
                        CL ( I ) = CL ( J )
                        CL ( J ) = ATEMP
                        ATEMP = CU ( I )
                        CU ( I ) = CU ( J )
                        CU ( J ) = ATEMP
                        GO TO 320
                     END IF
  310             CONTINUE
                  GO TO 330
               END IF
  320       CONTINUE
  330       CONTINUE
         END IF
      END IF

!  Initialize the performance counters and variables

 340  CONTINUE
      NC2OF = 0
      NC2OG = 0
      NC2OH = 0
      NC2CF = 0
      NC2CG = 0
      NC2CH = 0
      NHVPR = 0
      PNC = M
      STTIME = CPUTIM( DUM )
      SUTIME = STTIME - SUTIME

      RETURN

!  Non-executable statements.

 1000 FORMAT( I2, A8 )
 1001 FORMAT( 10I8 )
 1002 FORMAT( 2I8 )
 1010 FORMAT( ( 10I8 ) )
 1020 FORMAT( ( 1P, 4D16.8 ) )
 1030 FORMAT( ( 72L1 ) )
 1040 FORMAT( ( 8A10 ) )
 1080 FORMAT( 1P, 2D16.8 )
 1100 FORMAT( A8, 3I8 )
 1110 FORMAT( 1X, A6, /, ( 1X, 10I8 ) )
 1120 FORMAT( 1X, A6, /, ( 1X, 1P, 4D16.8 ) )
 1130 FORMAT( 1X, A6, /, ( 1X, 72L1 ) )
 1140 FORMAT( 1X, A6, /, ( 1X, 8A10 ) )
 1180 FORMAT( 1X, A6, /, 1P, 2D16.6 )
 2000 FORMAT( /, ' ** SUBROUTINE CSETUP: array length ', A6, &
              ' too small.', /, ' -- Miminimization abandoned.', &
              /, ' -- Increase the parameter ', A6, ' by at least ', I8, &
                 ' and restart.' )
 2010 FORMAT( /, ' ** SUBROUTINE CSETUP: ** Warning. The problem has', &
                 ' no general constraints. ', /, &
                 ' Other tools may be preferable' )
 2020 FORMAT( /, ' ** SUBROUTINE CSETUP: the problem uses no variables.' &
, ' Execution terminating ' )
 2030 FORMAT( /, ' ** SUBROUTINE CSETUP: the problem is vacuous.', &
                 ' Execution terminating ' )
 2040 FORMAT( ' ** SUBROUTINE CSETUP: Increase the size of IWK ' )

!  End of CSETUP.

      END
