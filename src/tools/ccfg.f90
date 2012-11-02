! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CCFG ( data, N, M, X, LC, C, JTRANS, LCJAC1, LCJAC2, CJAC, &
                         GRAD )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: N, M, LC, LCJAC1, LCJAC2
      REAL ( KIND = wp ) :: X( N ), C( LC ), CJAC( LCJAC1, LCJAC2 )
      LOGICAL :: JTRANS, GRAD

!  Compute the values of the constraint functions and their gradients
!  for constraints initially written in Standard Input Format (SIF).  
!  The Jacobian must be stored in a dense format.
! (Subroutine CSCFG performs the same calculations for a sparse Jacobian.)

!  CJAC  is a two-dimensional array of dimension ( LCJAC1, LCJAC2 )
!        which gives the value of the Jacobian matrix of the
!        constraint functions, or its transpose, evaluated at X.
!        If JTRANS is .TRUE., the i,j-th component of the array
!        will contain the i-th derivative of the j-th constraint
!        function. Otherwise, if JTRANS is .FALSE., the i,j-th
!        component of the array will contain the j-th derivative
!        of the i-th constraint function.

!  Based on the subroutines cfn.f and cgr.f by Nick Gould, which are
!  in turn based on the subroutine SBMIN by Conn, Gould and Toint.

!  Ingrid Bongartz
!  April 1992.


! ---------------------------------------------------------------------




! ---------------------------------------------------------------------



! ---------------------------------------------------------------------


! ---------------------------------------------------------------------

!  integer variables from the GLOBAL common block.


!  integer variables from the LOCAL common block.


!  integer variables from the PRFCTS common block.


!  local variables.

      INTEGER :: I, J, IEL, K, IG, II, IG1, L, LL, ICON
      INTEGER :: NIN, NVAREL, NELOW, NELUP, ISTRGV, IENDGV
      INTEGER :: ICNT, LLO, LLWRK, IFSTAT, IGSTAT
      LOGICAL :: NONTRV
!D    EXTERNAL           DSETVL, DSETVI, RANGE 
      REAL ( KIND = wp ) :: FTT, ONE, ZERO, GI, SCALEE
      PARAMETER ( ZERO = 0.0_wp, ONE = 1.0_wp )

      IF ( data%numcon == 0 ) RETURN

!  Check input parameters.

      IF( GRAD) THEN
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

!  Must identify which elements are included in constraints.
!  Use logical work vector to keep track of elements already included.
!  First ensure there is sufficient room in LOGI.

      LLO = GXEQX + data%ngng
      LLWRK = LLOGIC - LLO
      IF ( LLWRK < data%nelnum ) THEN
          IF ( IOUT > 0 ) WRITE( IOUT, 2020 ) data%nelnum - LLWRK 
          STOP
      END IF
      DO 410 I = 1, data%nelnum
         data%LOGIC( I ) = .FALSE.
  410 CONTINUE

!  Now identify elements in first M constraint groups.

      ICNT = 0
      DO 10 IG = 1, data%ng
         ICON = data%KNDOFC( IG )
         IF ( ICON > 0 .AND. ICON <= M ) THEN
            NELOW = data%ISTADG( IG )
            NELUP = data%ISTADG( IG + 1 ) - 1
            DO 20 II = NELOW, NELUP
               IEL = data%IELING( II )
               IF ( .NOT. data%LOGIC( IEL ) ) THEN
                  data%LOGIC( IEL ) = .TRUE.
                  ICNT = ICNT + 1
                  data%ICALCF( ICNT ) = IEL
               END IF
   20       CONTINUE
         END IF
   10 CONTINUE

!  Evaluate the element function values.

      CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), ICNT, &
                   data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                   data%IELVAR( 1 ), data%INTVAR( 1 ), &
                   data%ISTADH( 1 ), data%ISTEP( 1 ), &
                   data%ICALCF( 1 ),  &
                   data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                   data%lntvar, data%lintre, LFUVAL, data%lvscal, data%lepvlu,  &
                   1, IFSTAT )
      IF ( GRAD ) THEN

!  Evaluate the element function derivatives.

         CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), ICNT, &
                      data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                      data%IELVAR( 1 ), data%INTVAR( 1 ), &
                      data%ISTADH( 1 ), data%ISTEP( 1 ), &
                      data%ICALCF( 1 ),  &
                      data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                      data%lntvar, data%lintre, LFUVAL, data%lvscal, data%lepvlu,  &
                      2, IFSTAT )
      END IF

!  Compute the group argument values ft.

      DO 100 IG = 1, data%ng
         FTT = ZERO

!  Consider only those groups in the constraints.

         ICON = data%KNDOFC( IG )
         IF ( ICON > 0 .AND. ICON <= M ) THEN
            FTT = - data%B( IG )

!  Include the contribution from the linear element 
!  only if the variable belongs to the first N variables.

            DO 30 I = data%ISTADA( IG ), data%ISTADA( IG + 1 ) - 1
               J = data%ICNA( I )
               IF ( J <= N ) &
                  FTT = FTT + data%A( I ) * X( J )
   30       CONTINUE

!  Include the contributions from the nonlinear elements.

            DO 60 I = data%ISTADG( IG ), data%ISTADG( IG + 1 ) - 1
               FTT = FTT + &
                  data%ESCALE( I ) * data%FUVALS( data%IELING( I ) )
   60       CONTINUE

!  Record the derivatives of trivial groups.

            IF ( data%GXEQX( IG ) ) data%GVALS( data%ng + IG ) = ONE
         END IF
         data%FT( IG ) = FTT
  100 CONTINUE

!  Compute the group function values.

!  All group functions are trivial.

      IF ( data%altriv ) THEN
      CALL DCOPY( data%ng, data%FT( 1 ), 1, data%GVALS( 1 ), 1 )
      CALL DSETVL( data%ng, data%GVALS( data%ng + 1 ), 1, ONE )
      ELSE

!  Evaluate the group function values.
!  Evaluate groups belonging to the first M constraints only.

         ICNT = 0
         DO 400 IG = 1, data%ng
            ICON = data%KNDOFC( IG )
            IF ( ICON > 0 .AND. ICON <= M ) THEN
               ICNT = ICNT + 1
               data%ICALCF( ICNT ) = IG
            END IF 
  400    CONTINUE
         CALL GROUP ( data%GVALS( 1 ), data%ng, data%FT( 1 ), &
                      data%GPVALU( 1 ), ICNT, &
                      data%ITYPEG( 1 ), data%ISTGP( 1 ), &
                      data%ICALCF( 1 ),  &
                      data%lcalcg, data%ng1, data%lcalcg, data%lcalcg, data%lgpvlu, &
                      .FALSE., IGSTAT )
      END IF

!  Compute the constraint function values.

      DO 110 IG = 1, data%ng
         I = data%KNDOFC( IG )
         IF ( I > 0 .AND. I <= M ) THEN
            IF ( data%GXEQX( IG ) ) THEN
               C( I ) = data%GSCALE( IG ) * data%FT( IG )
            ELSE
               C( I ) = data%GSCALE( IG ) * data%GVALS( IG )
            END IF
         END IF
  110 CONTINUE
!  Increment the constraint function evaluation counter

      data%nc2cf = data%nc2cf + data%pnc
      IF ( GRAD ) THEN
!  Increment the constraint gradient evaluation counter

      data%nc2cg = data%nc2cg + data%pnc

!  Evaluate the group derivative values.

         IF ( .NOT. data%altriv ) CALL GROUP ( data%GVALS( 1 ), data%ng, &
               data%FT( 1 ), data%GPVALU( 1 ), ICNT, &
               data%ITYPEG( 1 ), data%ISTGP( 1 ), &
               data%ICALCF( 1 ), &
               data%lcalcg, data%ng1, data%lcalcg, data%lcalcg, data%lgpvlu, &
               .TRUE., IGSTAT )

!  Compute the gradient values.  Initialize the Jacobian as zero.

         DO 120 J = 1, N
            DO 115 I = 1, M
               IF ( JTRANS ) THEN
                  CJAC( J, I ) = ZERO
               ELSE
                  CJAC( I, J ) = ZERO
               END IF
  115       CONTINUE
  120    CONTINUE

!  Consider the IG-th group.

         DO 290 IG = 1, data%ng
            ICON = data%KNDOFC( IG )

!  Consider only those groups in the first M constraints.

            IF ( ICON == 0 .OR. ICON > M ) GO TO 290
            IG1 = IG + 1
            ISTRGV = data%IWORK( data%lstagv + IG )
            IENDGV = data%IWORK( data%lstagv + IG1 ) - 1
            NELOW = data%ISTADG( IG )
            NELUP = data%ISTADG( IG1 ) - 1
            NONTRV = .NOT. data%GXEQX( IG )

!  Compute the first derivative of the group.

            GI = data%GSCALE( IG )
            IF ( NONTRV ) GI = GI  * data%GVALS( data%ng + IG )

!  The group has nonlinear elements.

            IF ( NELOW <= NELUP ) THEN
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

!  Include contributions from the first N variables only.

                  IF ( LL <= N ) THEN
                     IF ( JTRANS ) THEN
                        CJAC( LL, ICON ) = GI * data%WRK( LL )
                     ELSE
                        CJAC( ICON, LL ) = GI * data%WRK( LL )
                     END IF
                  END IF
  190          CONTINUE

!  The group has only linear elements.

            ELSE

!  Allocate a gradient.

!DIR$ IVDEP
               DO 210 K = data%ISTADA( IG ), data%ISTADA( IG1 ) - 1
                  LL = data%ICNA( K )

!  Include contributions from the first N variables only.

                  IF ( LL <= N ) THEN
                     IF ( JTRANS ) THEN
!                            with only linear elements.
                        CJAC( LL, ICON ) = GI * data%A( K )
                     ELSE
!                            with only linear elements.
                        CJAC( ICON, LL ) = GI * data%A( K )
                     END IF
                  END IF
  210          CONTINUE
            END IF
  290    CONTINUE
      END IF
      RETURN

! Non-executable statements.

 2000 FORMAT( ' ** SUBROUTINE CCFG: Increase the leading dimension', &
              ' of CJAC ' ) 
 2010 FORMAT( ' ** SUBROUTINE CCFG: Increase the second dimension', &
              ' of CJAC ' )
 2020 FORMAT( /  ' ** SUBROUTINE CCFG: array length LLOGIC too small.' &
              /  ' -- Minimization abandoned.' &
              /  ' -- Increase the parameter LLOGIC by at least ', I8, &
                 ' and restart.' )

!  end of CCFG.

      END

