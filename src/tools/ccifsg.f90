! ( Last modified on 10 Sepc 2004 at 17:01:38 )
!  Correction: 10/Sep/2004: Incorrect argument to RANGES removed
      SUBROUTINE CCIFSG ( N, ICON, X, CI, NNZGCI, LGCI, GCI, &
                          INDVAR, GRAD )
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: N, ICON, NNZGCI, LGCI
      LOGICAL :: GRAD
      INTEGER :: INDVAR( LGCI )
      REAL ( KIND = wp ) :: CI
      REAL ( KIND = wp ) :: X( N ), GCI( LGCI )

!  Evaluate constraint function ICON and possibly its gradient,
!  for constraints initially written in Standard Input Format (SIF).
!  The constraint gradient is stored as a sparse vector in array GCI.
!  The j-th entry of GCI gives the value of the partial derivative
!  of constraint ICON with respect to variable INDVAR( j ).
!  The number of nonzeros in vector GCI is given by NNZGCI.
! (Subroutine CCIFG performs the same calculations for a dense
!  constraint gradient vector.)

!  Ingrid Bongartz
!  September 1994.

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

! ---------------------------------------------------------------------

!  End of parameters which might be changed by the user.

! ---------------------------------------------------------------------

      INTEGER :: IWK ( LIWK )
      REAL ( KIND = wp ) :: WK ( LWK )
      REAL ( KIND = wp ) :: FUVALS( LFUVAL )
      LOGICAL :: LOGI ( LLOGIC )
      CHARACTER ( LEN = 10 ) :: CHA ( LCHARA )

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
      INTEGER :: IOUT
      COMMON / OUTPUT /  IOUT
      INTEGER :: NUMVAR, NUMCON
      COMMON / DIMS /    NUMVAR, NUMCON
      SAVE             / GLOBAL /, / LOCAL /, / CHARA /, / OUTPUT /, &
                       / DIMS   /, / PRFCTS /

!  local variables.

      INTEGER :: I, J, IEL, K, IG, II, IG1, L, LL, NELING
      INTEGER :: NIN, NVAREL, NELOW, NELUP, ISTRGV, IENDGV
      INTEGER :: IFSTAT, IGSTAT
      LOGICAL :: NONTRV
!D    EXTERNAL           DSETVL, DSETVI, RANGE 
      REAL ( KIND = wp ) :: FTT, ONE, ZERO, GI, SCALEE
      PARAMETER ( ZERO = 0.0_wp, ONE = 1.0_wp )

!  Return if there are no constraints.

      IF ( NUMCON == 0 ) RETURN

!  Check input parameters.

      IF ( ICON <= 0 ) THEN
         WRITE( IOUT, 2000 )
         STOP
      END IF

!  Find group index IG of constraint ICON.

      IG = 0
      DO 70 I = 1, NG
         IF ( IWK( KNDOFC + I ) == ICON ) THEN
           IG = I
           GOTO 80
         END IF
   70 CONTINUE
   80 IF ( IG == 0 ) THEN
         WRITE( IOUT, 2000 )
         STOP
      END IF

!  Determine nonlinear elements in group IG.
!  Record their indices in IWK( ICALCF ).

      NELOW = IWK( ISTADG + IG )
      NELUP = IWK( ISTADG + IG + 1 ) - 1
      NELING = NELUP - NELOW + 1
      J = NELOW - 1
      DO 10 I = 1, NELING
         J = J + 1
         IWK( ICALCF + I ) = IWK( IELING + J )
   10 CONTINUE

!  Evaluate the element function values.

      CALL ELFUN ( FUVALS, X, WK( EPVALU + 1 ), NELING, &
                   IWK( ITYPEE + 1 ), IWK( ISTAEV + 1 ), &
                   IWK( IELVAR + 1 ), IWK( INTVAR + 1 ), &
                   IWK( ISTADH + 1 ), IWK( ISTEP + 1 ), &
                   IWK( ICALCF + 1 ),  &
                   LINTRE, LSTAEV, LELVAR, LNTVAR, LSTADH,  &
                   LNTVAR, LINTRE, LFUVAL, LVSCAL, LEPVLU,  &
                   1, IFSTAT )

!  Compute the group argument value FTT.
!  Consider only the group associated with constraint ICON.

      FTT = - WK( B + IG )

!  Include contributions from the linear element 
!  only if the variable belongs to the first N variables.

      DO 30 I = IWK( ISTADA + IG ), IWK( ISTADA + IG + 1 ) - 1
         J = IWK( ICNA + I )
         IF ( J <= N ) &
            FTT = FTT + WK( A + I ) * X( J )
   30 CONTINUE

!  Include the contributions from the nonlinear elements.

      DO 60 I = NELOW, NELUP
         FTT = FTT + &
                WK( ESCALE + I ) * FUVALS( IWK( IELING + I ) )
   60 CONTINUE
      WK( FT + IG ) = FTT

!  If IG is a trivial group, record the function value and derivative.

      IF ( LOGI( GXEQX + IG ) ) THEN
         WK( GVALS + IG ) = WK( FT + IG )
         WK( GVALS + NG + IG ) = ONE
      ELSE

!  Otherwise, evaluate group IG.

         CALL GROUP ( WK ( GVALS + 1 ), NG, WK( FT + 1 ), &
                      WK ( GPVALU + 1 ), 1, &
                      IWK( ITYPEG + 1 ), IWK( ISTGP + 1 ), &
                      IG, LCALCG, NG1, LCALCG, LCALCG, LGPVLU, &
                      .FALSE., IGSTAT )
      END IF

!  Compute the constraint function value.

      IF ( LOGI( GXEQX + IG ) ) THEN
         CI = WK( GSCALE + IG ) * WK( FT + IG )
      ELSE 
         CI = WK( GSCALE + IG ) * WK( GVALS + IG )
      END IF
      IF ( GRAD ) THEN

!  Evaluate the element function derivatives.

         CALL ELFUN ( FUVALS, X, WK( EPVALU + 1 ), NELING, &
                      IWK( ITYPEE + 1 ), IWK( ISTAEV + 1 ), &
                      IWK( IELVAR + 1 ), IWK( INTVAR + 1 ), &
                      IWK( ISTADH + 1 ), IWK( ISTEP + 1 ), &
                      IWK( ICALCF + 1 ),  &
                      LINTRE, LSTAEV, LELVAR, LNTVAR, LSTADH,  &
                      LNTVAR, LINTRE, LFUVAL, LVSCAL, LEPVLU,  &
                      2, IFSTAT )

!  Evaluate the group derivative values.

         IF ( .NOT. LOGI( GXEQX + IG ) ) &
            CALL GROUP ( WK ( GVALS + 1 ), NG, WK( FT + 1 ), &
                         WK ( GPVALU + 1 ), 1, &
                         IWK( ITYPEG + 1 ), IWK( ISTGP + 1 ), &
                         IG, LCALCG, NG1, LCALCG, LCALCG, LGPVLU, &
                         .TRUE., IGSTAT )

!  Compute the gradient values.  Initialize the gradient vector as zero.

         NNZGCI = 0 
         DO 120 J = 1, LGCI
            GCI( J ) = ZERO
  120    CONTINUE

!  Consider only group IG.

         IG1 = IG + 1
         ISTRGV = IWK( LSTAGV + IG )
         IENDGV = IWK( LSTAGV + IG1 ) - 1
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
  130             CONTINUE
               ELSE

!  The IEL-th element has no internal representation.

!DIR$ IVDEP
                  DO 140 I = 1, NVAREL
                     J = IWK( IELVAR + L )
                     WK( WRK + J ) = WK( WRK + J ) + &
                                        SCALEE * FUVALS( K )
                     K = K + 1
                     L = L + 1
  140             CONTINUE
               END IF
  150       CONTINUE

!  Include the contribution from the linear element.

!DIR$ IVDEP
            DO 160 K = IWK( ISTADA + IG ), &
                               IWK( ISTADA + IG1 ) - 1
               J = IWK( ICNA + K )
               WK( WRK + J ) = WK( WRK + J ) + WK( A + K )
  160       CONTINUE

!  Allocate a gradient.

!DIR$ IVDEP
            DO 190 I = ISTRGV, IENDGV
               LL = IWK( LSVGRP + I )

!  Include contributions from the first N variables only.

               IF ( LL <= N ) THEN
                  NNZGCI = NNZGCI + 1
                  GCI ( NNZGCI ) = GI * WK( WRK + LL )
                  INDVAR( NNZGCI ) = LL
               END IF
  190       CONTINUE

!  The group has only linear elements.

         ELSE

!  Include the contribution from the linear element.

!DIR$ IVDEP
            DO 210 K = IWK( ISTADA + IG ),IWK( ISTADA + IG1 ) - 1
               J = IWK( ICNA + K )
               WK( WRK + J ) = WK( WRK + J ) + WK( A + K )
  210       CONTINUE

!  Allocate a gradient.

!DIR$ IVDEP
            DO 220 I = ISTRGV, IENDGV
               LL = IWK( LSVGRP + I )

!  Include contributions from the first N variables only.

               IF ( LL <= N ) THEN
                  NNZGCI = NNZGCI + 1
                  GCI ( NNZGCI ) = GI * WK( WRK + LL )
                  INDVAR( NNZGCI ) = LL
               END IF
  220       CONTINUE
         END IF
      END IF

!  Update the counters for the report tool.

      NC2CF = NC2CF + 1
      IF ( GRAD ) NC2CG = NC2CG + 1
      RETURN

! Non-executable statements.

 2000 FORMAT( ' ** SUBROUTINE CCIFSG: invalid constraint index ICON ' )

!  end of CCIFSG.

      END
