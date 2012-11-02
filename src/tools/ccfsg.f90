! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CCFSG ( N, M, X, LC, C, NNZJ, LCJAC, CJAC, &
                         INDVAR, INDFUN, GRAD )
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: N, M, LC, NNZJ, LCJAC
      LOGICAL :: GRAD
      INTEGER :: INDVAR( LCJAC ), INDFUN( LCJAC )
      REAL ( KIND = wp ) :: X( N ), C( LC ), CJAC ( LCJAC )

!  Compute the values of the constraint functions and their gradients
!  for constraints initially written in Standard Input Format (SIF).
!  The Jacobian must be stored in a sparse format.
! (Subroutine CCFG performs the same calculations for a dense Jacobian.)

!  CJAC  is an array which gives the values of the nonzeros of the
!        general constraint functions evaluated at X and V.
!        The i-th entry of CJAC gives the value of the derivative
!        with respect to variable INDVAR(i) of constraint function 
!        INDFUN(i) (i.e., INDFUN(i) = j > 0 indicates the j-th
!        general constraint function).

!  Based on the subroutines cfn.f and csgr.f by Nick Gould, which are
!  in turn based on the subroutine SBMIN by Conn, Gould and Toint.

!  Ingrid Bongartz
!  April 1992.

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

!  integer variables from the PRFCTS common block.

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
      INTEGER :: NUMVAR, NUMCON
      COMMON / DIMS /    NUMVAR, NUMCON
      INTEGER :: IOUT
      COMMON / OUTPUT /  IOUT
      SAVE             / GLOBAL /, / LOCAL /, / CHARA /, / DIMS /, &
                       / PRFCTS /, / OUTPUT /

!  local variables.

      INTEGER :: I, J, IEL, K, IG, II, IG1, L, LL, ICON, ICNT
      INTEGER :: NIN, NVAREL, NELOW, NELUP, ISTRGV, IENDGV
      INTEGER :: LLO, LLWRK, IFSTAT, IGSTAT
      LOGICAL :: NONTRV
!D    EXTERNAL           DSETVL, DSETVI, RANGE 
      REAL ( KIND = wp ) :: FTT, ONE, ZERO, GI, SCALEE
      PARAMETER ( ZERO = 0.0_wp, ONE = 1.0_wp )

      IF ( NUMCON == 0 ) RETURN

!  Must identify which elements are included in constraints.
!  Use logical work vector to keep track of elements already included.
!  First ensure there is sufficient room in LOGI.

      LLO = GXEQX + NGNG
      LLWRK = LLOGIC - LLO
      IF ( LLWRK < NELNUM ) THEN
          IF ( IOUT > 0 ) WRITE( IOUT, 2010 ) NELNUM - LLWRK 
          STOP
      END IF
      DO 410 I = 1, NELNUM
         LOGI( LLO + I ) = .FALSE.
  410 CONTINUE

!  Now identify elements in first M constraint groups.

      ICNT = 0
      DO 10 IG = 1, NG
         ICON = IWK( KNDOFC + IG )
         IF ( ICON > 0 .AND. ICON <= M ) THEN
            NELOW = IWK( ISTADG + IG )
            NELUP = IWK( ISTADG + IG + 1 ) - 1
            DO 20 II = NELOW, NELUP
               IEL = IWK( IELING + II )
               IF ( .NOT. LOGI( LLO + IEL ) ) THEN
                  LOGI( LLO + IEL ) = .TRUE.
                  ICNT = ICNT + 1
                  IWK( ICALCF + ICNT ) = IEL
               END IF
   20       CONTINUE
         END IF
   10 CONTINUE

!  Evaluate the element function values.

      CALL ELFUN ( FUVALS, X, WK( EPVALU + 1 ), ICNT, &
                   IWK( ITYPEE + 1 ), IWK( ISTAEV + 1 ), &
                   IWK( IELVAR + 1 ), IWK( INTVAR + 1 ), &
                   IWK( ISTADH + 1 ), IWK( ISTEP + 1 ), &
                   IWK( ICALCF + 1 ),  &
                   LINTRE, LSTAEV, LELVAR, LNTVAR, LSTADH,  &
                   LNTVAR, LINTRE, LFUVAL, LVSCAL, LEPVLU,  &
                   1, IFSTAT )
      IF ( GRAD ) THEN

!  Evaluate the element function derivatives.

         CALL ELFUN ( FUVALS, X, WK( EPVALU + 1 ), ICNT, &
                      IWK( ITYPEE + 1 ), IWK( ISTAEV + 1 ), &
                      IWK( IELVAR + 1 ), IWK( INTVAR + 1 ), &
                      IWK( ISTADH + 1 ), IWK( ISTEP + 1 ), &
                      IWK( ICALCF + 1 ),  &
                      LINTRE, LSTAEV, LELVAR, LNTVAR, LSTADH,  &
                      LNTVAR, LINTRE, LFUVAL, LVSCAL, LEPVLU,  &
                      2, IFSTAT )
      END IF

!  Compute the group argument values ft.

      DO 100 IG = 1, NG
         FTT = ZERO

!  Consider only those groups in the constraints.

         ICON = IWK( KNDOFC + IG )
         IF ( ICON > 0 .AND. ICON <= M ) THEN
            FTT = - WK( B + IG )

!  Include contributions from the linear element 
!  only if the variable belongs to the first N variables.

            DO 30 I = IWK( ISTADA + IG ), IWK( ISTADA + IG + 1 ) - 1
               J = IWK( ICNA + I )
               IF ( J <= N ) &
                  FTT = FTT + WK( A + I ) * X( J )
   30       CONTINUE

!  Include the contributions from the nonlinear elements.

            DO 60 I = IWK( ISTADG + IG ), IWK( ISTADG + IG + 1 ) - 1
               FTT = FTT + &
                      WK( ESCALE + I ) * FUVALS( IWK( IELING + I ) )
   60       CONTINUE

!  Record derivatives of trivial groups.

            IF ( LOGI( GXEQX + IG ) ) WK( GVALS + NG + IG ) = ONE 
         END IF
         WK( FT + IG ) = FTT
  100 CONTINUE

!  Compute the group function values.

!  All group functions are trivial.

      IF ( ALTRIV ) THEN
      CALL DCOPY( NG, WK( FT + 1 ), 1, WK( GVALS + 1 ), 1 )
      CALL DSETVL( NG, WK( GVALS + NG + 1 ), 1, ONE )
      ELSE

!  Evaluate the group function values.
!  Evaluate groups belonging to the first M constraints only.

         ICNT = 0
         DO 400 IG = 1, NG
            ICON = IWK( KNDOFC + IG )
            IF ( ICON > 0 .AND. ICON <= M ) THEN
               ICNT = ICNT + 1
               IWK( ICALCF + ICNT ) = IG
            END IF 
  400    CONTINUE
         CALL GROUP ( WK ( GVALS + 1 ), NG, WK( FT + 1 ), &
                      WK ( GPVALU + 1 ), ICNT, &
                      IWK( ITYPEG + 1 ), IWK( ISTGP + 1 ), &
                      IWK( ICALCF + 1 ), &
                      LCALCG, NG1, LCALCG, LCALCG, LGPVLU, &
                      .FALSE., IGSTAT )
      END IF

!  Compute the constraint function values.

      DO 110 IG = 1, NG
         I = IWK( KNDOFC + IG )
         IF ( I > 0 .AND. I <= M ) THEN
            IF ( LOGI( GXEQX + IG ) ) THEN
               C( I ) = WK( GSCALE + IG ) * WK( FT + IG )
            ELSE 
               C( I ) = WK( GSCALE + IG ) * WK( GVALS + IG )
            END IF
         END IF
  110 CONTINUE
      IF ( GRAD ) THEN

!  Evaluate the group derivative values.

         IF ( .NOT. ALTRIV ) CALL GROUP ( WK ( GVALS + 1 ), NG, &
               WK( FT + 1 ), WK ( GPVALU + 1 ), ICNT, &
               IWK( ITYPEG + 1 ), IWK( ISTGP + 1 ), &
               IWK( ICALCF + 1 ), &
               LCALCG, NG1, LCALCG, LCALCG, LGPVLU, &
               .TRUE., IGSTAT )

!  Compute the gradient values.  Initialize the Jacobian as zero.

         NNZJ = 0 
         DO 120 J = 1, LCJAC
            CJAC( J ) = ZERO
  120    CONTINUE

!  Consider the IG-th group.

         DO 290 IG = 1, NG
            ICON = IWK( KNDOFC + IG )

!  Consider only those groups in the first M constraints.

            IF ( ICON == 0 .OR. ICON > M ) GO TO 290
            IG1 = IG + 1
            ISTRGV = IWK( LSTAGV + IG )
            IENDGV = IWK( LSTAGV + IG1 ) - 1
            NELOW = IWK( ISTADG + IG )
            NELUP = IWK( ISTADG + IG1 ) - 1
            NONTRV = .NOT. LOGI( GXEQX + IG )

!  Compute the first derivative of the group.

            GI = WK( GSCALE + IG )
            IF ( NONTRV ) GI = GI  * WK( GVALS + NG + IG )
      CALL DSETVI( IENDGV - ISTRGV + 1, WK( WRK + 1 ), &
                         IWK( LSVGRP + ISTRGV ), ZERO )

!  The group has nonlinear elements.

            IF ( NELOW <= NELUP ) THEN

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

!  Include contributions from the first N variables only.

                  IF ( LL <= N ) THEN
                     NNZJ = NNZJ + 1
                     IF ( NNZJ <= LCJAC ) THEN
                        CJAC ( NNZJ ) = GI * WK( WRK + LL )
                        INDFUN( NNZJ ) = ICON
                        INDVAR( NNZJ ) = LL
                     END IF
                  END IF
  190          CONTINUE

!  The group has only linear elements.

            ELSE
!                             linear element improved. 26 lines replace 19

!  Include the contribution from the linear element.

!DIR$ IVDEP
               DO 210 K = IWK( ISTADA + IG ),IWK( ISTADA + IG1 ) - 1
                  J = IWK( ICNA + K )
                  WK( WRK + J ) = WK( WRK + J ) + WK( A + K )
  210          CONTINUE

!  Allocate a gradient.

!DIR$ IVDEP
               DO 220 I = ISTRGV, IENDGV
                  LL = IWK( LSVGRP + I )

!  Include contributions from the first N variables only.

                  IF ( LL <= N ) THEN
                     NNZJ = NNZJ + 1
                     IF ( NNZJ <= LCJAC ) THEN
                        CJAC ( NNZJ ) = GI * WK( WRK + LL )
                        INDFUN( NNZJ ) = ICON
                        INDVAR( NNZJ ) = LL
                     END IF
                  END IF
  220          CONTINUE
            END IF
  290    CONTINUE

!  Verify that the Jacobian can fit in the allotted space

         IF ( NNZJ > LCJAC ) THEN
            IF ( IOUT > 0 ) WRITE( IOUT, 2000 ) NNZJ - LCJAC 
            STOP
         END IF
      END IF

!  Update the counters for the report tool.

      NC2CF = NC2CF + PNC
      IF ( GRAD ) NC2CG = NC2CG + PNC
      RETURN

!  Non-executable statements.

 2000 FORMAT( /  ' ** SUBROUTINE CCFSG: array length LCJAC too small.' &
              /  ' -- Minimization abandoned.', &
              /  ' -- Increase the parameter LCJAC by at least ', I8, &
                 ' and restart.' )
 2010 FORMAT( /  ' ** SUBROUTINE CCFSG: array length LLOGIC too small.' &
              /  ' -- Minimization abandoned.', &
              /  ' -- Increase the parameter LLOGIC by at least ', I8, &
                 ' and restart.' )

!  end of CSCFG.

      END

