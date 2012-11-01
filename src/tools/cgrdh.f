! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CGRDH ( N, M, X, GRLAGF, LV, V, G, &
                         JTRANS, LCJAC1, LCJAC2, CJAC, LH1, H )
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: N, M, LV, LH1, LCJAC1, LCJAC2
      LOGICAL :: GRLAGF, JTRANS
      REAL ( KIND = wp ) :: X ( N ), G ( N ), V ( LV ), &
                         CJAC ( LCJAC1, LCJAC2 ), &
                         H ( LH1, N )

!  Compute both the gradients of the objective, or Lagrangian, and
!  general constraint functions and the Hessian matrix of the
!  Lagrangian function of a problem initially written in Standard
!  Input Format (SIF).

!  G	 is an array which gives the value of the gradient of
!	 the objective function evaluated at X (GRLAGF = .FALSE.)
!        of of the Lagrangian function evaluated at X and V
! (GRLAGF = .TRUE.),

!  CJAC	 is a two-dimensional array of dimension ( LCJAC1, LCJAC2 )
!	 which gives the value of the Jacobian matrix of the
!	 constraint functions, or its transpose, evaluated at X.
!	 If JTRANS is .TRUE., the i,j-th component of the array
!        will contain the i-th derivative of the j-th constraint
!        function. Otherwise, if JTRANS is .FALSE., the i,j-th
!        component of the array will contain the j-th derivative
!        of the i-th constraint function.

!  H     is a two-dimensional array which gives the value of the
!        Hessian matrix of the Lagrangian function evaluated at
!        X and V. The i,j-th component of the array will contain the
!        derivative with respect to variables X(i) and X(j).

!  Based on the minimization subroutine LANCELOT/SBMIN
!  by Conn, Gould and Toint.

!  Nick Gould, for CGT productions,
!  November 1991.

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

!  integer variables from the GLOBAL common block.

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

!  integer variables from the LOCAL common block.

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
      COMMON / PRFCTS /  NC2OF, NC2OG, NC2OH, NC2CF, NC2CG, NC2CH, &
                         NHVPR, PNC, SUTIME, STTIME
      INTEGER :: IOUT
      COMMON / OUTPUT /  IOUT
      INTEGER :: NUMVAR, NUMCON
      COMMON / DIMS /    NUMVAR, NUMCON
      SAVE             / GLOBAL /, / LOCAL /, / CHARA /, / OUTPUT /, &
                       / DIMS   /, / PRFCTS /

!  Local variables

      INTEGER :: LH, LWKH, LIWKH, LIRNH, LJCNH, ICON
      INTEGER :: LNXTRW, LINXTR, INFORM, NNZH, IENDGV
      INTEGER :: I, J, IEL, K, IG, II, IG1, L, JJ, LL, LIH
      INTEGER :: NIN, NVAREL, NELOW, NELUP, ISTRGV
      INTEGER :: IFSTAT, IGSTAT
      LOGICAL :: NONTRV
!D    EXTERNAL           DSETVL, DSETVI, RANGE 
      REAL ( KIND = wp ) :: FTT, ONE, ZERO, GI, SCALEE, GII
      PARAMETER ( ZERO = 0.0_wp, ONE = 1.0_wp )

!  Check input parameters.

!                            dimension-checking.
      IF ( NUMCON > 0 ) THEN
         IF ( JTRANS ) THEN
            IF ( LCJAC1 < N .OR. LCJAC2 < M ) THEN
               IF ( LCJAC1 < N ) WRITE( IOUT, 2020 )
               IF ( LCJAC2 < M ) WRITE( IOUT, 2030 )
               STOP
            END IF
         ELSE
            IF ( LCJAC1 < M .OR. LCJAC2 < N ) THEN
               IF ( LCJAC1 < M ) WRITE( IOUT, 2020 )
               IF ( LCJAC2 < N ) WRITE( IOUT, 2030 )
               STOP
            END IF
         END IF
      END IF
      IF ( LH1 < N ) THEN
         WRITE( IOUT, 2040 )
         STOP
      END IF

!  there are non-trivial group functions.

      DO 10 I = 1, MAX( NELNUM, NG )
        IWK( ICALCF + I ) = I
   10 CONTINUE

!  evaluate the element function values.

      CALL ELFUN ( FUVALS, X, WK( EPVALU + 1 ), NELNUM, &
                   IWK( ITYPEE + 1 ), IWK( ISTAEV + 1 ), &
                   IWK( IELVAR + 1 ), IWK( INTVAR + 1 ), &
                   IWK( ISTADH + 1 ), IWK( ISTEP + 1 ), &
                   IWK( ICALCF + 1 ),  &
                   LINTRE, LSTAEV, LELVAR, LNTVAR, LSTADH,  &
                   LNTVAR, LINTRE, LFUVAL, LVSCAL, LEPVLU,  &
                   1, IFSTAT )

!  evaluate the element function values.

      CALL ELFUN ( FUVALS, X, WK( EPVALU + 1 ), NELNUM, &
                   IWK( ITYPEE + 1 ), IWK( ISTAEV + 1 ), &
                   IWK( IELVAR + 1 ), IWK( INTVAR + 1 ), &
                   IWK( ISTADH + 1 ), IWK( ISTEP + 1 ), &
                   IWK( ICALCF + 1 ),  &
                   LINTRE, LSTAEV, LELVAR, LNTVAR, LSTADH,  &
                   LNTVAR, LINTRE, LFUVAL, LVSCAL, LEPVLU,  &
                   3, IFSTAT )

!  compute the group argument values ft.

      DO 40 IG = 1, NG
         FTT = - WK( B + IG )

!  include the contribution from the linear element.

         DO 20 J = IWK( ISTADA + IG ), IWK( ISTADA + IG + 1 ) - 1
            FTT = FTT + WK( A + J ) * X( IWK( ICNA + J ) )
   20    CONTINUE

!  include the contributions from the nonlinear elements.

         DO 30 J = IWK( ISTADG + IG ), IWK( ISTADG + IG + 1 ) - 1
            FTT = FTT + WK( ESCALE + J ) * FUVALS( IWK( IELING + J ) )
   30    CONTINUE
         WK( FT + IG ) = FTT

!  Record the derivatives of trivial groups.

         IF ( LOGI( GXEQX + IG ) ) THEN
            WK( GVALS + NG + IG ) = ONE
            WK( GVALS + 2 * NG + IG ) = ZERO
         END IF
   40 CONTINUE

!  evaluate the group derivative values.

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

!  For unconstrained problems, skip construction of WK( LWKSTR + IG ) 
!  and skip specialized construction of gradient and Jacobian.  
!  Call ELGRD instead.

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
   80    CONTINUE

!  Compute the gradient values. Initialize the gradient and
!  Jacobian (or its transpose) as zero.

         DO 120 J = 1, N
            G( J ) = ZERO
            DO 110 I = 1, M
               IF ( JTRANS ) THEN
                  CJAC( J, I ) = ZERO
               ELSE
                  CJAC( I, J ) = ZERO
               END IF
  110       CONTINUE
  120    CONTINUE

!  Consider the IG-th group.

         DO 290 IG = 1, NG
            IG1 = IG + 1
            ICON = IWK( KNDOFC + IG )
            ISTRGV = IWK( LSTAGV + IG )
            IENDGV = IWK( LSTAGV + IG1 ) - 1
            NELOW = IWK( ISTADG + IG )
            NELUP = IWK( ISTADG + IG1 ) - 1
            NONTRV = .NOT. LOGI( GXEQX + IG )

!  Compute the first derivative of the group.

            GI = WK( GSCALE + IG )
            GII = WK( LWKSTR + IG )
            IF ( NONTRV ) THEN
               GI = GI  * WK( GVALS + NG + IG )
               GII = GII * WK( GVALS + NG + IG )
            END IF

!  This is the first gradient evaluation or the group has nonlinear
!  elements.

            IF ( FIRSTG .OR. NELOW <= NELUP ) THEN
      CALL DSETVI( IENDGV - ISTRGV + 1, WK( WRK + 1 ), &
                            IWK( LSVGRP + ISTRGV ), ZERO )

!  Loop over the group's nonlinear elements.

               DO 150 II = NELOW, NELUP
                  IEL = IWK( IELING + II )
                  K = IWK( INTVAR + IEL )
                  L = IWK( ISTAEV + IEL )
                  NVAREL = IWK( ISTAEV + IEL + 1 ) - L
                  SCALEE = WK( ESCALE + II )
                  IF ( LOGI( INTREP + IEL ) ) THEN

!  The IEL-th element has an internal representation.

                     NIN = IWK( INTVAR + IEL + 1 ) - K
                     CALL RANGE ( IEL, .TRUE., FUVALS( K ), &
                                  WK( WRK + N + 1 ), NVAREL, NIN, &
                                  IWK( ITYPEE + IEL ), &
                                  NIN, NVAREL )
!DIR$ IVDEP
                     DO 130 I = 1, NVAREL
                        J = IWK( IELVAR + L )
                        WK( WRK + J ) = WK( WRK + J ) + &
                                           SCALEE * WK( WRK + N + I )
                        L = L + 1
  130                CONTINUE
                  ELSE

!  The IEL-th element has no internal representation.

!DIR$ IVDEP
                     DO 140 I = 1, NVAREL
                        J = IWK( IELVAR + L )
                        WK( WRK + J ) = WK( WRK + J ) + &
                                           SCALEE * FUVALS( K )
                        K = K + 1
                        L = L + 1
  140                CONTINUE
                  END IF
  150          CONTINUE

!  Include the contribution from the linear element.

!DIR$ IVDEP
               DO 160 K = IWK( ISTADA + IG ), &
                                  IWK( ISTADA + IG1 ) - 1
                  J = IWK( ICNA + K )
                  WK( WRK + J ) = WK( WRK + J ) + WK( A + K )
  160          CONTINUE

!  Allocate a gradient.

!DIR$ IVDEP
               DO 190 I = ISTRGV, IENDGV
                  LL = IWK( LSVGRP + I )

!  The group belongs to the objective function.

                  IF ( ICON == 0 ) THEN
                     G( LL ) = G( LL ) + GI * WK( WRK + LL )

!  The group defines a constraint.

                  ELSE
                     IF ( JTRANS ) THEN
                        CJAC( LL, ICON ) = GI * WK( WRK + LL )
                     ELSE
                        CJAC( ICON, LL ) = GI * WK( WRK + LL )
                     END IF
                     IF ( GRLAGF ) &
                        G( LL ) = G( LL ) + GII * WK( WRK + LL )
                  END IF

!  If the group is non-trivial, also store the nonzero entries of the
!  gradient of the function in GRJAC.

                  IF ( NONTRV ) THEN
                     JJ = IWK( LSTAJC + LL )
                     FUVALS( LGRJAC + JJ ) = WK ( WRK + LL )

!  Increment the address for the next nonzero in the column of
!  the jacobian for variable LL.

                     IWK( LSTAJC + LL ) = JJ + 1
                  END IF
  190          CONTINUE

!  This is not the first gradient evaluation and there is only a linear
!  element.

            ELSE

!  Allocate a gradient.

!DIR$ IVDEP
               DO 210 K = IWK( ISTADA + IG ), IWK( ISTADA + IG1 ) - 1
                  LL = IWK( ICNA + K )

!  The group belongs to the objective function.

                  IF ( ICON == 0 ) THEN
                     G( LL ) = G( LL ) + GI * WK( A + K )

!  The group defines a constraint.

                  ELSE
                     IF ( JTRANS ) THEN
                        CJAC( LL, ICON ) =            GI * WK( A + K )
                     ELSE
                        CJAC( ICON, LL ) =            GI * WK( A + K )
                     END IF
                     IF ( GRLAGF ) G( LL ) = G( LL ) + GII * WK( A + K )
                  END IF
  210          CONTINUE

!  The group is non-trivial; increment the starting addresses for
!  the groups used by each variable in the (unchanged) linear
!  element to avoid resetting the nonzeros in the jacobian.

               IF ( NONTRV ) THEN
!DIR$ IVDEP
                  DO 220 I = ISTRGV, IENDGV
                     LL = IWK( LSVGRP + I )
                     IWK( LSTAJC + LL ) = IWK( LSTAJC + LL ) + 1
  220             CONTINUE
               END IF
            END IF
  290    CONTINUE

!  Reset the starting addresses for the lists of groups using
!  each variable to their values on entry.

         DO 300 I = N, 2, - 1
            IWK( LSTAJC + I ) = IWK( LSTAJC + I - 1 )
  300    CONTINUE
         IWK( LSTAJC + 1 ) = 1
      ELSE

!  Compute the gradient value.

      CALL DELGRD( N, NG, FIRSTG, IWK( ICNA + 1 ), LICNA, &
                      IWK( ISTADA + 1 ), LSTADA, IWK( IELING + 1 ), &
                      LELING, IWK( ISTADG + 1 ), LSTADG, &
                      IWK( ITYPEE + 1 ), LINTRE, &
                      IWK( ISTAEV + 1 ), LSTAEV, IWK( IELVAR + 1 ), &
                      LELVAR, IWK( INTVAR + 1 ), LNTVAR, &
                      IWK( LSVGRP + 1 ), &
                      LNVGRP, IWK( LSTAJC + 1 ), LNSTJC, &
                      IWK( LSTAGV + 1 ), LNSTGV, WK( A + 1 ), LA, &
                      WK( GVALS + NG + 1 ), LGVALS, &
                      FUVALS, LNGUVL, FUVALS( LGGFX + 1 ), &
                      WK( GSCALE + 1 ), LGSCAL, &
                      WK( ESCALE + 1 ), LESCAL, FUVALS( LGRJAC + 1 ), &
                      LNGRJC, WK( WRK + 1 ), WK( WRK + N + 1 ), MAXSEL, &
                      LOGI( GXEQX + 1 ), LGXEQX, &
                      LOGI( INTREP + 1 ), LINTRE, RANGE )

!  Store the gradient value.

         DO 400 I = 1, N
            G( I ) = FUVALS( LGGFX + I )
  400    CONTINUE
      END IF
      FIRSTG = .FALSE.

!  Define the real work space needed for ASMBL.
!  Ensure that there is sufficient space.

!                            for unconstrained problems.
      IF ( NUMCON > 0 ) THEN
         LWKH = LWK2 - N - 3 * MAXSEL - NG
      ELSE
         LWKH = LWK2 - N - 3 * MAXSEL
      END IF
      IF ( LWKH <= 0 ) THEN
         IF ( IOUT > 0 ) WRITE( IOUT, 2000 )
         STOP
      END IF

!  Define the integer work space needed for ASMBL.
!  Ensure that there is sufficient space.

      LIWKH = LIWK2 - N
      LH = MIN( LWKH, ( LIWKH - 3 * N ) / 4 )
      LINXTR = LH + N
      IF ( LH <= 0 ) THEN
         IF ( IOUT > 0 ) WRITE( IOUT, 2010 )
         STOP
      END IF

!  Set starting addresses for partitions of the integer workspace.

      LIH = LH
      LIRNH = 0
      LJCNH = LIRNH + LIH
      LNXTRW = LJCNH + LIH
      DO 310 I = 1, N
         IWK( IVAR + I ) = I
  310 CONTINUE

!  Assemble the Hessian.

      IF ( NUMCON > 0 ) THEN
      CALL DASMBL( N, NG, MAXSEL, N, LH, LIH, NNZH, &
                   N, IWK( IVAR + 1), IWK( ISTADH + 1 ), LSTADH, &
                   IWK( ICNA + 1 ), LICNA, &
                   IWK( ISTADA + 1 ), LSTADA, IWK( INTVAR + 1 ), LNTVAR, &
                   IWK( IELVAR + 1 ), LELVAR, IWK( IELING + 1 ), LELING, &
                   IWK( ISTADG + 1 ), LSTADG, IWK( ISTAEV + 1 ), LSTAEV, &
                   IWK( LSTAGV + 1 ), LNSTGV, IWK( LSVGRP + 1 ), LNVGRP, &
                   IWK( LSEND + LIRNH + 1 ), IWK ( LSEND + LJCNH + 1 ), &
                   IWK( LSEND + LNXTRW + 1 ), LINXTR, &
                   IWK( LSEND + LIWKH + 1 ), N, &
                   WK( A + 1 ), LA, FUVALS, LNGUVL, FUVALS, LNHUVL, &
                   WK( GVALS + NG + 1 ), WK( GVALS + 2 * NG + 1 ), &
                   WK( LWKSTR + 1 ), WK( ESCALE + 1 ), LESCAL, &
                   WK( LWKSTR + NG + 1 ), WK( LWKSTR + LWKH + NG + 1 ), &
                   LWK2 - LWKH, LOGI( GXEQX + 1 ), &
                   LGXEQX, LOGI( INTREP + 1 ), LINTRE, &
                   IWK( ITYPEE + 1 ), LINTRE, &
                   RANGE, 1, IOUT, .FALSE., I, INFORM, &
                   .FALSE., .FALSE. )
      ELSE
      CALL DASMBL( N, NG, MAXSEL, N, LH, LIH, NNZH, &
                   N, IWK( IVAR + 1), IWK( ISTADH + 1 ), LSTADH, &
                   IWK( ICNA + 1 ), LICNA, &
                   IWK( ISTADA + 1 ), LSTADA, IWK( INTVAR + 1 ), LNTVAR, &
                   IWK( IELVAR + 1 ), LELVAR, IWK( IELING + 1 ), LELING, &
                   IWK( ISTADG + 1 ), LSTADG, IWK( ISTAEV + 1 ), LSTAEV, &
                   IWK( LSTAGV + 1 ), LNSTGV, IWK( LSVGRP + 1 ), LNVGRP, &
                   IWK( LSEND + LIRNH + 1 ), IWK ( LSEND + LJCNH + 1 ), &
                   IWK( LSEND + LNXTRW + 1 ), LINXTR, &
                   IWK( LSEND + LIWKH + 1 ), N, &
                   WK( A + 1 ), LA, FUVALS, LNGUVL, FUVALS, LNHUVL, &
                   WK( GVALS + NG + 1 ), WK( GVALS + 2 * NG + 1 ), &
                   WK( GSCALE + 1 ), WK( ESCALE + 1 ), LESCAL, &
                   WK( LWKSTR + 1 ), WK( LWKSTR + LWKH + 1 ), &
                   LWK2 - LWKH, LOGI( GXEQX + 1 ), &
                   LGXEQX, LOGI( INTREP + 1 ), LINTRE, &
                   IWK( ITYPEE + 1 ), LINTRE, &
                   RANGE, 1, IOUT, .FALSE., I, INFORM, &
                   .FALSE., .FALSE. )
      END IF

!  Check that there is sufficient integer workspace.

      IF ( INFORM > 0 ) THEN
         IF ( IOUT > 0 ) WRITE( IOUT, 2010 )
         STOP
      END IF

!  Initialize the dense Hessian matrix.

      DO 420 J = 1, N
         DO 410 I = 1, N
            H( I, J ) = ZERO
  410    CONTINUE
  420 CONTINUE

!  Transfer the matrix from co-ordinate to dense storage and
!  symmetrize the martix.

      IF ( NUMCON > 0 ) THEN
         DO 430 K = 1, NNZH
            I = IWK( LSEND + LIRNH + K )
            J = IWK( LSEND + LJCNH + K )
            H( I, J ) = WK( LWKSTR + NG + K )
            H( J, I ) = WK( LWKSTR + NG + K )
  430    CONTINUE
      ELSE
         DO 440 K = 1, NNZH
            I = IWK( LSEND + LIRNH + K )
            J = IWK( LSEND + LJCNH + K )
            H( I, J ) = WK( LWKSTR + K )
            H( J, I ) = WK( LWKSTR + K )
  440    CONTINUE
      END IF

!  Update the counters for the report tool.

      NC2OG = NC2OG + 1
      NC2CG = NC2CG + PNC
      NC2OH = NC2OH + 1
      NC2CH = NC2CH + PNC
      RETURN

! Non-executable statements.

 2000 FORMAT( ' ** SUBROUTINE CGRDH: Increase the size of WK ' )
 2010 FORMAT( ' ** SUBROUTINE CGRDH: Increase the size of IWK ' )
 2020 FORMAT( ' ** SUBROUTINE CGRDH: Increase the leading dimension', &
              ' of CJAC ' )
 2030 FORMAT( ' ** SUBROUTINE CGRDH: Increase the second dimension', &
              ' of CJAC ' )
 2040 FORMAT( ' ** SUBROUTINE CGRDH: Increase the leading dimension', &
              ' of H ' )

!  end of CGRDH.

      END
