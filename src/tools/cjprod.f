! ( Last modified on 10 Sepc 2004 at 16:55:38 )
!  Correction: 10/Sep/2004: undeclared integers variables declared
      SUBROUTINE CJPROD( N, M, GOTJ, JTRANS, X, V, LV, R, LR )
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: N, M, LV, LR
      LOGICAL :: GOTJ, JTRANS
      REAL ( KIND = wp ) :: X( N ), V( LV ), R( LR )

!  Compute the matrix-vector product between the Jacobian matrix
!  of the constraints (JTRANS = .FALSE.), or its transpose 
! (JTRANS = .TRUE.) for the problem, and a given vector P. 
!  The result is placed in R. If GOTJ is .TRUE. the first derivatives 
!  are assumed to have already been computed. If the user is unsure, 
!  set GOTJ = .FALSE. the first time a product is required with the 
!  Jacobian evaluated at X. X is not used if GOTJ = .TRUE.

!  Nick Gould, for GOT/CUTEr productions.
!  June, 2003.

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

      INTEGER :: I, IG, J, ICON, K, IG1, II
      INTEGER :: L, IEL, NVAREL, NIN
      INTEGER :: IFSTAT, IGSTAT
      REAL ( KIND = wp ) :: ZERO, ONE, FTT, PROD, SCALEE
      PARAMETER ( ZERO = 0.0_wp, ONE = 1.0_wp )
      EXTERNAL :: RANGE 
      IF ( NUMCON == 0 ) RETURN

!  Check input data.

      IF ( ( JTRANS .AND. LV < M ) .OR.  &
 ( .NOT. JTRANS .AND. LV < N ) ) THEN
         IF ( IOUT > 0 ) WRITE( IOUT, 2010 )
         STOP
      END IF
      IF ( ( JTRANS .AND. LR < N ) .OR.  &
 ( .NOT. JTRANS .AND. LR < M ) ) THEN
         IF ( IOUT > 0 ) WRITE( IOUT, 2020 )
         STOP
      END IF

!  There are non-trivial group functions.

      IF ( .NOT. GOTJ ) THEN
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
            END IF
   70    CONTINUE

!  Evaluate the group derivative values.

         IF ( .NOT. ALTRIV ) CALL GROUP ( WK ( GVALS + 1 ), NG, &
               WK( FT + 1 ), WK ( GPVALU + 1 ), NG, &
               IWK( ITYPEG + 1 ), IWK( ISTGP + 1 ), &
               IWK( ICALCF + 1 ), &
               LCALCG, NG1, LCALCG, LCALCG, LGPVLU, &
               .TRUE., IGSTAT )
      END IF

!  Ensure that there is sufficient space.

      IF ( LWK2 < N ) THEN
         IF ( IOUT > 0 ) WRITE( IOUT, 2000 )
         STOP
      END IF

!  Form the product r = J(transpose) v

      IF ( JTRANS ) THEN

!  Initialize R

         DO 110 I = 1, N
            R( I ) = ZERO
  110    CONTINUE

!  Consider the IG-th group.

         DO 190 IG = 1, NG
            ICON = IWK( KNDOFC + IG )
            IF ( ICON > 0 ) THEN
               IG1 = IG + 1


!  Compute the product of v(i) with the (scaled) group derivative

               IF ( LOGI( GXEQX + IG ) ) THEN
                  PROD = V( ICON ) * WK( GSCALE + IG )
               ELSE
                  PROD = V( ICON ) * WK( GSCALE + IG ) * WK( GVALS + NG + IG )
               END IF

!  Loop over the group's nonlinear elements.

               DO 150 II = IWK( ISTADG + IG ),  &
                           IWK( ISTADG + IG1 ) - 1
                  IEL = IWK( IELING + II )
                  K = IWK( INTVAR + IEL )
                  L = IWK( ISTAEV + IEL )
                  NVAREL = IWK( ISTAEV + IEL + 1 ) - L
                  SCALEE = WK( ESCALE + II ) * PROD
                  IF ( LOGI( INTREP + IEL ) ) THEN

!  The IEL-th element has an internal representation.

                     NIN = IWK( INTVAR + IEL + 1 ) - K
                     CALL RANGE ( IEL, .TRUE., FUVALS( K ), &
                                  WK( WRK + 1 ), NVAREL, NIN, &
                                  IWK( ITYPEE + IEL ), &
                                  NIN, NVAREL )
!DIR$ IVDEP
                     DO 130 I = 1, NVAREL
                        J = IWK( IELVAR + L )
                        R( J ) = R( J ) + SCALEE * WK( WRK + I )
                        L = L + 1
  130                CONTINUE
                  ELSE

!  The IEL-th element has no internal representation.

!DIR$ IVDEP
                     DO 140 I = 1, NVAREL
                        J = IWK( IELVAR + L )
                        R( J ) = R( J ) + SCALEE * FUVALS( K )
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
                  R( J ) = R( J ) + WK( A + K ) * PROD
  160          CONTINUE
            END IF
  190    CONTINUE

!  Form the product r = J v

      ELSE

!  Consider the IG-th group.

         DO 290 IG = 1, NG
            ICON = IWK( KNDOFC + IG )
            IF ( ICON > 0 ) THEN
               IG1 = IG + 1
               PROD = ZERO

!  Compute the first derivative of the group.

!  Loop over the group's nonlinear elements.

               DO 250 II = IWK( ISTADG + IG ),  &
                           IWK( ISTADG + IG1 ) - 1
                  IEL = IWK( IELING + II )
                  K = IWK( INTVAR + IEL )
                  L = IWK( ISTAEV + IEL )
                  NVAREL = IWK( ISTAEV + IEL + 1 ) - L
                  SCALEE = WK( ESCALE + II )
                  IF ( LOGI( INTREP + IEL ) ) THEN

!  The IEL-th element has an internal representation.

                     NIN = IWK( INTVAR + IEL + 1 ) - K
                     CALL RANGE ( IEL, .TRUE., FUVALS( K ), &
                                  WK( WRK + 1 ), NVAREL, NIN, &
                                  IWK( ITYPEE + IEL ), &
                                  NIN, NVAREL )
!DIR$ IVDEP
                     DO 230 I = 1, NVAREL
                        PROD = PROD + V( IWK( IELVAR + L ) ) * &
                                       SCALEE * WK( WRK + I )
                        L = L + 1
  230                CONTINUE
                  ELSE

!  The IEL-th element has no internal representation.

!DIR$ IVDEP
                     DO 240 I = 1, NVAREL
                        PROD = PROD + V( IWK( IELVAR + L ) ) * &
                                   SCALEE * FUVALS( K )
                        K = K + 1
                        L = L + 1
  240                CONTINUE
                  END IF
  250          CONTINUE

!  Include the contribution from the linear element.

!DIR$ IVDEP
               DO 260 K = IWK( ISTADA + IG ), &
                          IWK( ISTADA + IG1 ) - 1
                  PROD = PROD + V( IWK( ICNA + K ) ) * WK( A + K )
  260          CONTINUE

!  Multiply the product by the (scaled) group derivative

               IF ( LOGI( GXEQX + IG ) ) THEN
                  R( ICON ) = PROD * WK( GSCALE + IG )
               ELSE
                  R( ICON ) = PROD * WK( GSCALE + IG ) * WK( GVALS + NG + IG )
               END IF
            END IF
  290    CONTINUE
      END IF

!  Update the counters for the report tool.

      IF ( .NOT. GOTJ ) THEN
         NC2OG = NC2OG + 1
         NC2CG = NC2CG + PNC
      END IF
      RETURN

! Non-executable statements.

 2000 FORMAT( ' ** SUBROUTINE CJPROD: Increase the size of WK ' )
 2010 FORMAT( ' ** SUBROUTINE CJPROD: Increase the size of V ' )
 2020 FORMAT( ' ** SUBROUTINE CJPROD: Increase the size of R ' )

!  end of CJPROD.

      END
