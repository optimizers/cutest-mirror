! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CGR ( data, N, M, X, GRLAGF, LV, V, &
                         G, JTRANS, LCJAC1, LCJAC2, CJAC )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: N, M, LV, LCJAC1, LCJAC2
      LOGICAL :: GRLAGF, JTRANS
      REAL ( KIND = wp ) :: X( N ), G( N ), V( LV ), CJAC( LCJAC1, LCJAC2 )

!  Compute both the gradients of the objective, or Lagrangian, and
!  general constraint functions of a problem initially written in
!  Standard Input Format (SIF).

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

!  Based on the minimization subroutine data%laNCELOT/SBMIN
!  by Conn, Gould and Toint.

!  Nick Gould, for CGT productions,
!  November 1991.


! ---------------------------------------------------------------------




! ---------------------------------------------------------------------



! ---------------------------------------------------------------------


! ---------------------------------------------------------------------

!  integer variables from the GLOBAL common block.


!  integer variables from the LOCAL common block.


!  Integer variables from the PRFCTS common block.


!  local variables.

      INTEGER :: I, J, IEL, K, IG, II, IG1, L, JJ, LL, ICON
      INTEGER :: NIN, NVAREL, NELOW, NELUP, ISTRGV, IENDGV
      INTEGER :: IFSTAT, IGSTAT
      LOGICAL :: NONTRV
!D    EXTERNAL           DSETVL, DSETVI, RANGE 
      REAL ( KIND = wp ) :: FTT, ONE, ZERO, GI, SCALEE, GII
      PARAMETER ( ZERO = 0.0_wp, ONE = 1.0_wp )

!  Check input parameters.

!                            dimension-checking.
      IF ( data%numcon > 0 ) THEN
         IF ( JTRANS ) THEN
            IF ( LCJAC1 < N .OR. LCJAC2 < M ) THEN
               IF ( LCJAC1 < N ) WRITE( IOUT, 2000 )
               IF ( LCJAC2 < M ) WRITE( IOUT, 2010 )
               STOP
            END IF
         ELSE
            IF ( LCJAC1 < M .OR. LCJAC2 < N ) THEN
               IF ( LCJAC1 < M ) WRITE( IOUT, 2000 )
               IF ( LCJAC2 < N ) WRITE( IOUT, 2010 )
               STOP
            END IF
         END IF
      END IF

!  there are non-trivial group functions.

      DO 10 I = 1, MAX( data%nelnum, data%ng )
        data%ICALCF( I ) = I
   10 CONTINUE

!  evaluate the element function values.

      CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), data%nelnum, &
                   data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                   data%IELVAR( 1 ), data%INTVAR( 1 ), &
                   data%ISTADH( 1 ), data%ISTEP( 1 ), &
                   data%ICALCF( 1 ),  &
                   data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                   data%lntvar, data%lintre, LFUVAL, data%lvscal, data%lepvlu,  &
                   1, IFSTAT )

!  evaluate the element function values.

      CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), data%nelnum, &
                   data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                   data%IELVAR( 1 ), data%INTVAR( 1 ), &
                   data%ISTADH( 1 ), data%ISTEP( 1 ), &
                   data%ICALCF( 1 ),  &
                   data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                   data%lntvar, data%lintre, LFUVAL, data%lvscal, data%lepvlu,  &
                   2, IFSTAT )

!  compute the group argument values ft.

      DO 40 IG = 1, data%ng
         FTT = - data%B( IG )

!  include the contribution from the linear element.

         DO 20 J = data%ISTADA( IG ), data%ISTADA( IG + 1 ) - 1
            FTT = FTT + data%A( J ) * X( data%ICNA( J ) )
   20    CONTINUE

!  include the contributions from the nonlinear elements.

         DO 30 J = data%ISTADG( IG ), data%ISTADG( IG + 1 ) - 1
            FTT = FTT + data%ESCALE( J ) * data%FUVALS( data%IELING( J ) )
   30    CONTINUE
         data%FT( IG ) = FTT

!  Record the derivatives of trivial groups.

         IF ( data%GXEQX( IG ) ) data%GVALS( data%ng + IG ) = ONE
   40 CONTINUE

!  evaluate the group derivative values.

      IF ( .NOT. data%altriv ) CALL GROUP ( data%GVALS( 1 ), data%ng, &
            data%FT( 1 ), data%GPVALU( 1 ), data%ng, &
            data%ITYPEG( 1 ), data%ISTGP( 1 ), &
            data%ICALCF( 1 ), &
            data%lcalcg, data%ng1, data%lcalcg, data%lcalcg, data%lgpvlu, &
            .TRUE., IGSTAT )

!  For unconstrained problems, skip construction of gradient 
!  and Jacobian. Call ELGRD instead.

      IF ( data%numcon > 0 ) THEN

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

         DO 290 IG = 1, data%ng
            IG1 = IG + 1
            ICON = data%KNDOFC( IG )
            ISTRGV = data%IWORK( data%lstagv + IG )
            IENDGV = data%IWORK( data%lstagv + IG1 ) - 1
            NELOW = data%ISTADG( IG )
            NELUP = data%ISTADG( IG1 ) - 1
            NONTRV = .NOT. data%GXEQX( IG )

!  Compute the first derivative of the group.

            GI = data%GSCALE( IG )
            IF ( ICON == 0 ) THEN
               GII = GI
            ELSE
               IF ( GRLAGF ) GII = GI * V( data%KNDOFC( IG ) )
            END IF
            IF ( NONTRV ) THEN
               GI = GI  * data%GVALS( data%ng + IG )
               IF ( GRLAGF ) GII = GII * data%GVALS( data%ng + IG )
            END IF

!  This is the first gradient evaluation or the group has nonlinear
!  elements.

            IF ( data%firstg .OR. NELOW <= NELUP ) THEN
      CALL DSETVI( IENDGV - ISTRGV + 1, data%WRK( 1 ), &
                            data%IWORK( data%lsvgrp + ISTRGV ), ZERO )

!  Loop over the group's nonlinear elements.

               DO 150 II = NELOW, NELUP
                  IEL = data%IELING( II )
                  K = data%INTVAR( IEL )
                  L = data%ISTAEV( IEL )
                  NVAREL = data%ISTAEV( IEL + 1 ) - L
                  SCALEE = data%ESCALE( II )
                  IF ( data%INTREP( IEL ) ) THEN

!  The IEL-th element has an internal representation.

                     NIN = data%INTVAR( IEL + 1 ) - K
                     CALL RANGE ( IEL, .TRUE., data%FUVALS( K ), &
                                  data%WRK( N + 1 ), NVAREL, NIN, &
                                  data%ITYPEE( IEL ), &
                                  NIN, NVAREL )
!DIR$ IVDEP
                     DO 130 I = 1, NVAREL
                        J = data%IELVAR( L )
                        data%WRK( J ) = data%WRK( J ) + &
                                           SCALEE * data%WRK( N + I )
                        L = L + 1
  130                CONTINUE
                  ELSE

!  The IEL-th element has no internal representation.

!DIR$ IVDEP
                     DO 140 I = 1, NVAREL
                        J = data%IELVAR( L )
                        data%WRK( J ) = data%WRK( J ) + &
                                           SCALEE * data%FUVALS( K )
                        K = K + 1
                        L = L + 1
  140                CONTINUE
                  END IF
  150          CONTINUE

!  Include the contribution from the linear element.

!DIR$ IVDEP
               DO 160 K = data%ISTADA( IG ), &
                                  data%ISTADA( IG1 ) - 1
                  J = data%ICNA( K )
                  data%WRK( J ) = data%WRK( J ) + data%A( K )
  160          CONTINUE

!  Allocate a gradient.

!DIR$ IVDEP
               DO 190 I = ISTRGV, IENDGV
                  LL = data%IWORK( data%lsvgrp + I )

!  The group belongs to the objective function.

                  IF ( ICON == 0 ) THEN
                     G( LL ) = G( LL ) + GI * data%WRK( LL )

!  The group defines a constraint.

                  ELSE
                     IF ( JTRANS ) THEN
                        CJAC( LL, ICON ) =   GI * data%WRK( LL )
                     ELSE
                        CJAC( ICON, LL ) =   GI * data%WRK( LL )
                     END IF
                     IF ( GRLAGF ) &
                        G( LL ) = G( LL ) + GII * data%WRK( LL )
                  END IF

!  If the group is non-trivial, also store the nonzero entries of the
!  gradient of the function in GRJAC.

                  IF ( NONTRV ) THEN
                     JJ = data%IWORK( data%lstajc + LL )
                     data%FUVALS( data%lgrjac + JJ ) = data%WRK( LL )

!  Increment the address for the next nonzero in the column of
!  the jacobian for variable LL.

                     data%IWORK( data%lstajc + LL ) = JJ + 1
                  END IF
  190          CONTINUE

!  This is not the first gradient evaluation and there is only a linear
!  element.

            ELSE

!  Allocate a gradient.

!DIR$ IVDEP
               DO 210 K = data%ISTADA( IG ), data%ISTADA( IG1 ) - 1
                  LL = data%ICNA( K )

!  The group belongs to the objective function.

                  IF ( ICON == 0 ) THEN
                     G( LL ) = G( LL ) + GI * data%A( K )

!  The group defines a constraint.

                  ELSE
                     IF ( JTRANS ) THEN
                        CJAC( LL, ICON ) =            GI * data%A( K )
                     ELSE
                        CJAC( ICON, LL ) =            GI * data%A( K )
                     END IF
                     IF ( GRLAGF ) G( LL ) = G( LL ) + GII * data%A( K )
                  END IF
  210          CONTINUE

!  The group is non-trivial; increment the starting addresses for
!  the groups used by each variable in the (unchanged) linear
!  element to avoid resetting the nonzeros in the jacobian.

               IF ( NONTRV ) THEN
!DIR$ IVDEP
                  DO 220 I = ISTRGV, IENDGV
                     LL = data%IWORK( data%lsvgrp + I )
                     data%IWORK( data%lstajc + LL ) = data%IWORK( data%lstajc + LL ) + 1
  220             CONTINUE
               END IF
            END IF
  290    CONTINUE

!  Reset the starting addresses for the lists of groups using
!  each variable to their values on entry.

         DO 300 I = N, 2, - 1
            data%IWORK( data%lstajc + I ) = data%IWORK( data%lstajc + I - 1 )
  300    CONTINUE
         data%IWORK( data%lstajc + 1 ) = 1
      ELSE

!  Compute the gradient value.

      CALL DELGRD( N, data%ng, data%firstg, data%ICNA( 1 ), data%licna, &
                      data%ISTADA( 1 ), data%lstada, data%IELING( 1 ), &
                      data%leling, data%ISTADG( 1 ), data%lstadg, &
                      data%ITYPEE( 1 ), data%lintre, &
                      data%ISTAEV( 1 ), data%lstaev, data%IELVAR( 1 ), &
                      data%lelvar, data%INTVAR( 1 ), data%lntvar, &
                      data%IWORK( data%lsvgrp + 1 ), &
                      data%lnvgrp, data%IWORK( data%lstajc + 1 ), data%lnstjc, &
                      data%IWORK( data%lstagv + 1 ), data%lnstgv, data%A( 1 ), data%la, &
                      data%GVALS( data%ng + 1 ), data%lgvals, &
                      data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ), &
                      data%GSCALE( 1 ), data%lgscal, &
                      data%ESCALE( 1 ), data%lescal, data%FUVALS( data%lgrjac + 1 ), &
                      data%lngrjc, data%WRK( 1 ), data%WRK( N + 1 ), data%maxsel, &
                      data%GXEQX( 1 ), data%lgxeqx, &
                      data%INTREP( 1 ), data%lintre, RANGE )

!  Store the gradient value.

         DO 400 I = 1, N
            G( I ) = data%FUVALS( data%lggfx + I )
  400    CONTINUE
      END IF
      data%firstg = .FALSE.

!  Update the counters for the report tool.

      data%nc2og = data%nc2og + 1
      data%nc2cg = data%nc2cg + data%pnc
      RETURN

! Non-executable statements.

 2000 FORMAT( ' ** SUBROUTINE CGR: Increase the leading dimension', &
              ' of CJAC ' )
 2010 FORMAT( ' ** SUBROUTINE CGR: Increase the second dimension', &
              ' of CJAC ' )

!  end of CGR.

      END
