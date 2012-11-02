! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE COFG ( data, N, X, F, G, GRAD )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: N
      REAL ( KIND = wp ) :: F
      REAL ( KIND = wp ) :: X( N ), G( * )
      LOGICAL :: GRAD

!  Compute the value of the objective function and its gradient
!  for a function initially written in Standard Input Format (SIF).

!  G     is an array which gives the value of the gradient of
!        the objective function evaluated at X.
!        G(i) gives the partial derivative of the objective
!        function with respect to variable X(i).

!  Based on the subroutines cfn.f and cgr.f by Nick Gould, which are
!  in turn based on the subroutine SBMIN by Conn, Gould and Toint.

!  Ingrid Bongartz 
!  April 1992.


! ---------------------------------------------------------------------




! ---------------------------------------------------------------------



! ---------------------------------------------------------------------


! ---------------------------------------------------------------------

!  integer variables from the GLOBAL common block.

      INTEGER :: ESCALE, GSCALE, VSCALE, GVALS, XT, DGRAD 

!  integer variables from the LOCAL common block.


!  Integer variables from the PRFCTS common block.


!  local variables.

      INTEGER :: I, J, IEL, K, IG, II, IG1, L, LL, ICON
      INTEGER :: NIN, NVAREL, NELOW, NELUP, ISTRGV, IENDGV
      INTEGER :: ICNT, LLO, LLWRK, IFSTAT, IGSTAT
      LOGICAL :: NONTRV
!D    EXTERNAL           DSETVL, DSETVI, RANGE 
      REAL ( KIND = wp ) :: FTT, ONE, ZERO, GI, SCALEE
      PARAMETER ( ZERO = 0.0_wp, ONE = 1.0_wp ) 

!  Must identify which elements are included in objective function.
!  Use logical work vector to keep track of elements already included.
!  First ensure there is sufficient room in LOGI.

      LLO = GXEQX + data%ngng
      LLWRK = LLOGIC - LLO
      IF ( LLWRK < data%nelnum ) THEN
          IF ( IOUT > 0 ) WRITE( IOUT, 2000 ) data%nelnum - LLWRK 
          STOP
      END IF
      DO 410 I = 1, data%nelnum
         data%LOGIC( I ) = .FALSE.
  410 CONTINUE

!  Now identify elements in objective function groups.

      ICNT = 0
      DO 10 IG = 1, data%ng
         IF ( data%KNDOFC( IG ) == 0 ) THEN
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

!                            evaluation block if there are no constraints.
      IF ( data%numcon > 0 ) THEN
         DO 100 IG = 1, data%ng
            FTT = ZERO

!  Consider only those groups in the objective function. 

            IF ( data%KNDOFC( IG ) == 0 ) THEN
               FTT = - data%B( IG )

!  Include the contribution from the linear element 
!  only if the variable belongs to the first N variables.

               DO 30 I = data%ISTADA( IG ), data%ISTADA( IG + 1 ) - 1
                  J = data%ICNA( I ) 
                  IF ( J <= N )  &
                     FTT = FTT + data%A( I ) * X( J )
   30          CONTINUE

!  Include the contributions from the nonlinear elements.

               DO 60 I = data%ISTADG( IG ), data%ISTADG( IG + 1 ) - 1
                  FTT = FTT + &
                         data%ESCALE( I ) * data%FUVALS( data%IELING( I ) )
   60          CONTINUE

!  Record the derivatives of trivial groups.

               IF ( data%GXEQX( IG ) ) data%GVALS( data%ng + IG ) = ONE
            END IF
            data%FT( IG ) = FTT
  100    CONTINUE
      ELSE

!  There are no constraints, so we need not check data%KNDOFC( IG ).

         DO 300 IG = 1, data%ng
            FTT = - data%B( IG )

!  Include the contribution from the linear element
!  only if the variable belongs to the first N variables.

            DO 330 I = data%ISTADA( IG ), data%ISTADA( IG + 1 ) - 1
               J = data%ICNA( I )
               IF ( J <= N ) &
                  FTT = FTT + data%A( I ) * X( J )
  330       CONTINUE

!  Include the contributions from the nonlinear elements.

            DO 360 I = data%ISTADG( IG ), data%ISTADG( IG + 1 ) - 1
               FTT = FTT + &
                      data%ESCALE( I ) * data%FUVALS( data%IELING( I ) )
  360       CONTINUE

!  Record the derivatives of trivial groups.

            IF ( data%GXEQX( IG ) ) data%GVALS( data%ng + IG ) = ONE
            data%FT( IG ) = FTT
  300    CONTINUE
      END IF

!  Compute the group function values.

!  All group functions are trivial.

      IF ( data%altriv ) THEN
      CALL DCOPY( data%ng, data%FT( 1 ), 1, data%GVALS( 1 ), 1 )
      CALL DSETVL( data%ng, data%GVALS( data%ng + 1 ), 1, ONE )
      ELSE

!  Evaluate the group function values.
!  Evaluate groups belonging to the objective function only.

         ICNT = 0
         DO 400 IG = 1, data%ng
            IF ( data%KNDOFC( IG ) == 0 ) THEN
               ICNT = ICNT + 1
               data%ICALCF( ICNT ) = IG
            END IF 
  400    CONTINUE
         CALL GROUP ( data%GVALS( 1 ), data%ng, data%FT( 1 ), &
                      data%GPVALU( 1 ), ICNT, &
                      data%ITYPEG( 1 ), data%ISTGP( 1 ), &
                      data%ICALCF( 1 ), &
                      data%lcalcg, data%ng1, data%lcalcg, data%lcalcg, data%lgpvlu, &
                      .FALSE., IGSTAT )
      END IF

!  Compute the objective function value.

!                            block if there are no constraints.
      F = ZERO
      IF ( data%numcon > 0 ) THEN
         DO 110 IG = 1, data%ng
            IF ( data%KNDOFC( IG ) == 0 ) THEN
               IF ( data%GXEQX( IG ) ) THEN
                  F = F + data%GSCALE( IG ) * data%FT( IG )
               ELSE
                  F = F + data%GSCALE( IG ) * data%GVALS( IG )
               END IF
            END IF 
  110    CONTINUE
      ELSE

!  There are no constraints, so we need not check data%KNDOFC( IG ).

         DO 310 IG = 1, data%ng
            IF ( data%GXEQX( IG ) ) THEN
               F = F + data%GSCALE( IG ) * data%FT( IG )
            ELSE
               F = F + data%GSCALE( IG ) * data%GVALS( IG )
            END IF
  310    CONTINUE
      END IF
      IF ( GRAD ) THEN

!  Evaluate the group derivative values.

         IF ( .NOT. data%altriv ) CALL GROUP ( data%GVALS( 1 ), data%ng, &
               data%FT( 1 ), data%GPVALU( 1 ), ICNT, &
               data%ITYPEG( 1 ), data%ISTGP( 1 ), &
               data%ICALCF( 1 ), &
               data%lcalcg, data%ng1, data%lcalcg, data%lcalcg, data%lgpvlu, &
               .TRUE., IGSTAT )

!  Compute the gradient values. Initialize the gradient as zero.

         DO 120 J = 1, N
            G( J ) = ZERO
  120    CONTINUE   

!  Consider the IG-th group.

         DO 290 IG = 1, data%ng
            ICON = data%KNDOFC( IG )

!  Consider only those groups in the objective function. 

            IF ( ICON > 0 ) GO TO 290
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
      CALL DSETVI( IENDGV - ISTRGV + 1, data%WRK( 1 ),  &
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

!  Include the contributions from only the first N variables.

                  IF ( LL <= N ) &
                     G( LL ) = G( LL ) + GI * data%WRK( LL )
  190          CONTINUE

!  The group has only linear elements.

            ELSE

!  Allocate a gradient.

!DIR$ IVDEP
               DO 210 K = data%ISTADA( IG ), data%ISTADA( IG1 ) - 1
                  LL = data%ICNA( K )

!  Include the contributions from linear elements for only the first
!  N variables.

                  IF ( LL <= N ) &
                     G( LL ) = G( LL ) + GI * data%A( K )
  210          CONTINUE
            END IF
  290    CONTINUE
      ENDIF

!  Update the counters for the report tool.

      data%nc2of = data%nc2of + 1
      IF ( GRAD ) data%nc2og = data%nc2og + 1

      RETURN

!  Non-executable statements.

 2000 FORMAT( /  ' ** SUBROUTINE COFG: array length LLOGIC too small.' &
              /  ' -- Minimization abandoned.' &
              /  ' -- Increase the parameter LLOGIC by at least ', I8, &
                 ' and restart.' )

!  end of COFG.

      END

