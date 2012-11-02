! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CPROD( N, M, GOTH, X, LV, V, P, RESULT )
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: N, M, LV
      LOGICAL :: GOTH
      REAL ( KIND = wp ) :: X( N ), V ( LV ), P( N ), RESULT( N )

!  Compute the matrix-vector product between the Hessian matrix
!  of the Lagrangian function for the problem and  a given vector P.
!  The result is placed in RESULT. If GOTH is .TRUE. the second
!  derivatives are assumed to have already been computed. If
!  the user is unsure, set GOTH = .FALSE. the first time a product
!  is required with the Hessian evaluated at X and V. X and V are
!  not used if GOTH = .TRUE.

!  Based on the minimization subroutine LANCELOT/SBMIN
!  by Conn, Gould and Toint.

!  Nick Gould, for CGT productions.
!  November, 1991.

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
      INTEGER :: LWK2, MAXSIN, NINVAR, MAXSEL
      INTEGER :: NTYPE, NSETS, LSTYPE, LSSWTR, LSSIWT, LSIWTR
      INTEGER :: LSWTRA, LNTYPE, LNSWTR, LNSIWT, LNIWTR
      INTEGER :: LNWTRA, LSISET, LSSVSE, LNISET, LNSVSE
      LOGICAL :: ALTRIV, FIRSTG

!  Integer variables from the PRFCTS common block.

      INTEGER :: NC2OF, NC2OG, NC2OH, NC2CF, NC2CG, NC2CH
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
      INTEGER :: IOUT
      COMMON / OUTPUT /  IOUT
      COMMON / PRFCTS /  NC2OF, NC2OG, NC2OH, NC2CF, NC2CG, NC2CH, &
                         NHVPR, PNC, SUTIME, STTIME
      INTEGER :: NUMVAR, NUMCON
      COMMON / DIMS /    NUMVAR, NUMCON
      SAVE             / GLOBAL /, / LOCAL /, / CHARA /, / OUTPUT /, &
                       / DIMS   /, / PRFCTS /

!  Local variables

      INTEGER :: I, IG, J, NN, NBPROD, NNONNZ
      INTEGER :: LNWK, LNWKB, LNWKC, LWKB, LWKC
      INTEGER :: IFSTAT, IGSTAT
      REAL ( KIND = wp ) :: ZERO, ONE, FTT
      PARAMETER ( ZERO = 0.0_wp, ONE = 1.0_wp )
      EXTERNAL :: RANGE 

!  There are non-trivial group functions.

      IF ( .NOT. GOTH ) THEN
         DO 10 I = 1, MAX( NELNUM, NG )
           IWK( ICALCF + I ) = I
   10    CONTINUE

!  Evaluate the element function values.

         CALL ELFUN ( FUVALS, X, WK( EPVALU + 1 ), NELNUM, &
                      IWK( ITYPEE + 1 ), IWK( ISTAEV + 1 ), &
                      IWK( IELVAR + 1 ), IWK( INTVAR + 1 ), &
                      IWK( ISTADH + 1 ), IWK( ISTEP + 1 ), &
                      IWK( ICALCF + 1 ),  &
                      LINTRE, LSTAEV, LELVAR, LNTVAR, LSTADH,  &
                      LNTVAR, LINTRE, LFUVAL, LVSCAL, LEPVLU,  &
                      1, IFSTAT )

!  Evaluate the element function values.

         CALL ELFUN ( FUVALS, X, WK( EPVALU + 1 ), NELNUM, &
                      IWK( ITYPEE + 1 ), IWK( ISTAEV + 1 ), &
                      IWK( IELVAR + 1 ), IWK( INTVAR + 1 ), &
                      IWK( ISTADH + 1 ), IWK( ISTEP + 1 ), &
                      IWK( ICALCF + 1 ),  &
                      LINTRE, LSTAEV, LELVAR, LNTVAR, LSTADH,  &
                      LNTVAR, LINTRE, LFUVAL, LVSCAL, LEPVLU,  &
                      3, IFSTAT )

!  Compute the group argument values ft.

         DO 70 IG = 1, NG
            FTT = - WK( B + IG )

!  Include the contribution from the linear element.

            DO 30 J = IWK( ISTADA + IG ), IWK( ISTADA + IG + 1 ) - 1
               FTT = FTT + WK( A + J ) * X( IWK( ICNA + J ) )
   30       CONTINUE

!  Include the contributions from the nonlinear elements.

            DO 60 J = IWK( ISTADG + IG ), IWK( ISTADG + IG + 1 ) - 1
               FTT = FTT + WK( ESCALE + J ) * FUVALS( IWK( IELING + J))
   60       CONTINUE
            WK( FT + IG ) = FTT

!  Record the derivatives of trivial groups.

            IF ( LOGI( GXEQX + IG ) ) THEN
               WK( GVALS + NG + IG ) = ONE
               WK( GVALS + 2 * NG + IG ) = ZERO
            END IF
   70    CONTINUE

!  Evaluate the group derivative values.

         IF ( .NOT. ALTRIV ) CALL GROUP ( WK ( GVALS + 1 ), NG, &
               WK( FT + 1 ), WK ( GPVALU + 1 ), NG, &
               IWK( ITYPEG + 1 ), IWK( ISTGP + 1 ), &
               IWK( ICALCF + 1 ), &
               LCALCG, NG1, LCALCG, LCALCG, LGPVLU, &
               .TRUE., IGSTAT )

!  Define the real work space needed for ELGRD.
!  Ensure that there is sufficient space.

         IF ( LWK2 < NG ) THEN
            IF ( IOUT > 0 ) WRITE( IOUT, 2000 )
            STOP
         END IF
!                            for unconstrained problems.
         IF ( NUMCON > 0 ) THEN

!  Change the group weightings to include the contributions from
!  the Lagrange multipliers.

            DO 80 IG = 1, NG
               I = IWK( KNDOFC + IG )
               IF ( I == 0 ) THEN
                  WK( LWKSTR + IG ) = WK( GSCALE + IG )
               ELSE
                  WK( LWKSTR + IG ) = WK( GSCALE + IG ) * V( I )
               END IF
   80       CONTINUE

!  Compute the gradient value.

      CALL DELGRD( N, NG, FIRSTG, IWK( ICNA + 1 ), LICNA, &
                         IWK( ISTADA + 1 ), LSTADA, IWK( IELING + 1 ), &
                         LELING, IWK( ISTADG + 1 ), LSTADG, &
                         IWK( ITYPEE + 1 ), LINTRE, &
                         IWK( ISTAEV + 1 ), LSTAEV, IWK( IELVAR + 1 ), &
                         LELVAR, IWK( INTVAR + 1 ), LNTVAR, &
                         IWK( LSVGRP + 1 ), LNVGRP, IWK( LSTAJC + 1 ), &
                         LNSTJC, IWK( LSTAGV + 1 ), LNSTGV, &
                         WK( A + 1 ), LA, WK( GVALS + NG + 1 ), LGVALS, &
                         FUVALS, LNGUVL, FUVALS( LGGFX + 1 ), &
                         WK( LWKSTR + 1 ), NG, &
                         WK( ESCALE + 1 ), LESCAL, FUVALS( LGRJAC + 1 ), &
                         LNGRJC, WK( WRK + 1 ), WK( WRK + N + 1 ), &
                         MAXSEL, LOGI( GXEQX + 1 ), LGXEQX, &
                         LOGI( INTREP + 1 ), LINTRE, RANGE )
         ELSE

!  Compute the gradient value.

      CALL DELGRD( N, NG, FIRSTG, IWK( ICNA + 1 ), LICNA, &
                         IWK( ISTADA + 1 ), LSTADA, IWK( IELING + 1 ), &
                         LELING, IWK( ISTADG + 1 ), LSTADG, &
                         IWK( ITYPEE + 1 ), LINTRE, &
                         IWK( ISTAEV + 1 ), LSTAEV, IWK( IELVAR + 1 ), &
                         LELVAR, IWK( INTVAR + 1 ), LNTVAR, &
                         IWK( LSVGRP + 1 ), LNVGRP, IWK( LSTAJC + 1 ), &
                         LNSTJC, IWK( LSTAGV + 1 ), LNSTGV, &
                         WK( A + 1 ), LA, WK( GVALS + NG + 1 ), LGVALS, &
                         FUVALS, LNGUVL, FUVALS( LGGFX + 1 ), &
                         WK( GSCALE + 1 ), LGSCAL, &
                         WK( ESCALE + 1 ), LESCAL, FUVALS( LGRJAC + 1 ), &
                         LNGRJC, WK( WRK + 1 ), WK( WRK + N + 1 ), &
                         MAXSEL, LOGI( GXEQX + 1 ), LGXEQX, &
                         LOGI( INTREP + 1 ), LINTRE, RANGE )
         END IF
         FIRSTG = .FALSE.
      END IF

!  Ensure that the product involves all components of P.

      DO 100 I = 1, N
         IWK( IVAR + I ) = I
         IWK( LNNONZ + I ) = I
  100 CONTINUE

!  Initialize RESULT as the zero vector.

      CALL DSETVL( N, RESULT, 1, ZERO )

!  Define the real work space needed for HSPRD.
!  Ensure that there is sufficient space.

      NN = NINVAR + N
      LNWK = MAX( NG, MAXSEL )
      LNWKB = MAXSIN
      LNWKC = MAXSIN
      LWKB = LNWK
      LWKC = LWKB + LNWKB

!  Evaluate the product.

      IF ( NUMCON > 0 ) THEN
      CALL DHSPRD( N, NN, NG, NGEL, N, 1, N, NBPROD, NELNUM == 0, &
                   IWK( IVAR + 1 ), IWK( ISTAEV + 1 ), LSTAEV, &
                   IWK( ISTADH + 1 ), LSTADH, IWK( INTVAR + 1 ), &
                   LNTVAR, IWK( IELING + 1 ), LELING, IWK( IELVAR + 1 ), &
                   LELVAR, IWK( LSTAJC + 1 ), LNSTJC, IWK( LSELTS + 1 ), &
                   LNELTS, IWK( LSPTRS + 1 ), LNPTRS, IWK( LGCOLJ + 1 ),  &
                   LNGCLJ, IWK( LSLGRP + 1 ), LNLGRP, IWK( LSWKSP + 1 ),  &
                   LNWKSP, IWK( LSVGRP + 1 ), LNVGRP, IWK( LSTAGV + 1 ), &
                   LNSTGV, IWK( LVALJR + 1 ), LNVLJR, IWK( ITYPEE + 1 ),  &
                   LINTRE, NNONNZ, IWK( LNNONZ + 1 ), LNNNON, &
                   IWK( LIUSED + 1 ), LNIUSE, IWK( LNONZ2 + 1 ), &
                   LNNNO2, IWK( LSYMMH + 1 ), MAXSIN, P, RESULT, &
                   WK( GVALS + NG + 1 ), WK( GVALS + 2 * NG + 1 ), &
                   FUVALS( LGRJAC + 1 ), LNGRJC, WK( LWKSTR + 1 ), &
                   WK( ESCALE + 1 ), LESCAL, FUVALS, LNHUVL, &
                   WK( WRK + 1 ), LNWK, WK( WRK + LWKB + 1 ), &
                   LNWKB, WK( WRK + LWKC + 1 ), LNWKC, &
                   LOGI( GXEQX + 1 ), LGXEQX, LOGI( INTREP + 1 ), &
                   LINTRE, .TRUE., RANGE )
      ELSE
      CALL DHSPRD( N, NN, NG, NGEL, N, 1, N, NBPROD, NELNUM == 0, &
                   IWK( IVAR + 1 ), IWK( ISTAEV + 1 ), LSTAEV, &
                   IWK( ISTADH + 1 ), LSTADH, IWK( INTVAR + 1 ), &
                   LNTVAR, IWK( IELING + 1 ), LELING, IWK( IELVAR + 1 ), &
                   LELVAR, IWK( LSTAJC + 1 ), LNSTJC, IWK( LSELTS + 1 ), &
                   LNELTS, IWK( LSPTRS + 1 ), LNPTRS, IWK( LGCOLJ + 1 ),  &
                   LNGCLJ, IWK( LSLGRP + 1 ), LNLGRP, IWK( LSWKSP + 1 ),  &
                   LNWKSP, IWK( LSVGRP + 1 ), LNVGRP, IWK( LSTAGV + 1 ), &
                   LNSTGV, IWK( LVALJR + 1 ), LNVLJR, IWK( ITYPEE + 1 ),  &
                   LINTRE, NNONNZ, IWK( LNNONZ + 1 ), LNNNON, &
                   IWK( LIUSED + 1 ), LNIUSE, IWK( LNONZ2 + 1 ), &
                   LNNNO2, IWK( LSYMMH + 1 ), MAXSIN, P, RESULT, &
                   WK( GVALS + NG + 1 ), WK( GVALS + 2 * NG + 1 ), &
                   FUVALS( LGRJAC + 1 ), LNGRJC, WK( GSCALE + 1 ), &
                   WK( ESCALE + 1 ), LESCAL, FUVALS, LNHUVL, &
                   WK( WRK + 1 ), LNWK, WK( WRK + LWKB + 1 ), &
                   LNWKB, WK( WRK + LWKC + 1 ), LNWKC, &
                   LOGI( GXEQX + 1 ), LGXEQX, LOGI( INTREP + 1 ), &
                   LINTRE, .TRUE., RANGE )
      END IF

!  Update the counters for the report tool.

      NHVPR = NHVPR + 1
      IF ( .NOT. GOTH ) THEN
         NC2OH = NC2OH + 1
         NC2CH = NC2CH + PNC
      END IF
      RETURN

! Non-executable statements.

 2000 FORMAT( ' ** SUBROUTINE CPROD: Increase the size of WK ' )

!  end of CPROD.

      END


