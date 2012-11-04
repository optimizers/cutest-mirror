! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE COFG ( data, n, X, f, G, GRAD )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n
      REAL ( KIND = wp ) :: f
      REAL ( KIND = wp ) :: X( n ), G( * )
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

      INTEGER :: ESCALE, GSCALE, VSCALE, gvals, XT, DGRAD 

!  integer variables from the LOCAL common block.


!  Integer variables from the PRFCTS common block.


!  local variables.

      INTEGER :: i, j, iel, k, ig, ii, ig1, l, ll, icon
      INTEGER :: nin, nvarel, nelow, nelup, istrgv, iendgv
      INTEGER :: icnt, llo, llwrk, ifstat, igstat
      LOGICAL :: NONTRV
!D    EXTERNAL           SETVL, SETVI, RANGE 
      REAL ( KIND = wp ) :: ftt, one, zero, gi, scalee
      PARAMETER ( zero = 0.0_wp, one = 1.0_wp ) 

!  Must identify which elements are included in objective function.
!  Use logical work vector to keep track of elements already included.
!  First ensure there is sufficient room in LOGI.

      llo = gxeqx + data%ngng
      llwrk = llogic - llo
      IF ( llwrk < data%nelnum ) THEN
          IF ( iout > 0 ) WRITE( iout, 2000 ) data%nelnum - llwrk 
          STOP
      END IF
      DO 410 i = 1, data%nelnum
         data%LOGIC( i ) = .FALSE.
  410 CONTINUE

!  Now identify elements in objective function groups.

      icnt = 0
      DO 10 ig = 1, data%ng
         IF ( data%KNDOFC( ig ) == 0 ) THEN
            nelow = data%ISTADG( ig )
            nelup = data%ISTADG( ig + 1 ) - 1
            DO 20 ii = nelow, nelup
               iel = data%IELING( ii )
               IF ( .NOT. data%LOGIC( iel ) ) THEN
                  data%LOGIC( iel ) = .TRUE.
                  icnt = icnt + 1
                  data%ICALCF( icnt ) = iel
               END IF
   20       CONTINUE
         END IF
   10 CONTINUE

!  Evaluate the element function values.

      CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), icnt, &
                   data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                   data%IELVAR( 1 ), data%INTVAR( 1 ), &
                   data%ISTADH( 1 ), data%ISTEP( 1 ), &
                   data%ICALCF( 1 ),  &
                   data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                   data%lntvar, data%lintre, lfuval, data%lvscal, data%lepvlu,  &
                   1, ifstat )
      IF ( GRAD ) THEN

!  Evaluate the element function derivatives.

         CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), icnt, &
                      data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                      data%IELVAR( 1 ), data%INTVAR( 1 ), &
                      data%ISTADH( 1 ), data%ISTEP( 1 ), &
                      data%ICALCF( 1 ),  &
                   data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                   data%lntvar, data%lintre, lfuval, data%lvscal, data%lepvlu,  &
                   2, ifstat )
      END IF

!  Compute the group argument values ft.

!                            evaluation block if there are no constraints.
      IF ( data%numcon > 0 ) THEN
         DO 100 ig = 1, data%ng
            ftt = zero

!  Consider only those groups in the objective function. 

            IF ( data%KNDOFC( ig ) == 0 ) THEN
               ftt = - data%B( ig )

!  Include the contribution from the linear element 
!  only if the variable belongs to the first n variables.

               DO 30 i = data%ISTADA( ig ), data%ISTADA( ig + 1 ) - 1
                  j = data%ICNA( i ) 
                  IF ( j <= n )  &
                     ftt = ftt + data%A( i ) * X( j )
   30          CONTINUE

!  Include the contributions from the nonlinear elements.

               DO 60 i = data%ISTADG( ig ), data%ISTADG( ig + 1 ) - 1
                  ftt = ftt + &
                         data%ESCALE( i ) * data%FUVALS( data%IELING( i ) )
   60          CONTINUE

!  Record the derivatives of trivial groups.

               IF ( data%GXEQX( ig ) ) data%GVALS( data%ng + ig ) = one
            END IF
            data%FT( ig ) = ftt
  100    CONTINUE
      ELSE

!  There are no constraints, so we need not check data%KNDOFC( ig ).

         DO 300 ig = 1, data%ng
            ftt = - data%B( ig )

!  Include the contribution from the linear element
!  only if the variable belongs to the first n variables.

            DO 330 i = data%ISTADA( ig ), data%ISTADA( ig + 1 ) - 1
               j = data%ICNA( i )
               IF ( j <= n ) &
                  ftt = ftt + data%A( i ) * X( j )
  330       CONTINUE

!  Include the contributions from the nonlinear elements.

            DO 360 i = data%ISTADG( ig ), data%ISTADG( ig + 1 ) - 1
               ftt = ftt + &
                      data%ESCALE( i ) * data%FUVALS( data%IELING( i ) )
  360       CONTINUE

!  Record the derivatives of trivial groups.

            IF ( data%GXEQX( ig ) ) data%GVALS( data%ng + ig ) = one
            data%FT( ig ) = ftt
  300    CONTINUE
      END IF

!  Compute the group function values.

!  All group functions are trivial.

      IF ( data%altriv ) THEN
      CALL DCOPY( data%ng, data%FT( 1 ), 1, data%GVALS( 1 ), 1 )
      CALL SETVL( data%ng, data%GVALS( data%ng + 1 ), 1, one )
      ELSE

!  Evaluate the group function values.
!  Evaluate groups belonging to the objective function only.

         icnt = 0
         DO 400 ig = 1, data%ng
            IF ( data%KNDOFC( ig ) == 0 ) THEN
               icnt = icnt + 1
               data%ICALCF( icnt ) = ig
            END IF 
  400    CONTINUE
         CALL GROUP ( data%GVALS( 1 ), data%ng, data%FT( 1 ), &
                      data%GPVALU( 1 ), icnt, &
                      data%ITYPEG( 1 ), data%ISTGP( 1 ), &
                      data%ICALCF( 1 ), &
                      data%lcalcg, data%ng1, data%lcalcg, data%lcalcg, data%lgpvlu, &
                      .FALSE., igstat )
      END IF

!  Compute the objective function value.

!                            block if there are no constraints.
      f = zero
      IF ( data%numcon > 0 ) THEN
         DO 110 ig = 1, data%ng
            IF ( data%KNDOFC( ig ) == 0 ) THEN
               IF ( data%GXEQX( ig ) ) THEN
                  f = f + data%GSCALE( ig ) * data%FT( ig )
               ELSE
                  f = f + data%GSCALE( ig ) * data%GVALS( ig )
               END IF
            END IF 
  110    CONTINUE
      ELSE

!  There are no constraints, so we need not check data%KNDOFC( ig ).

         DO 310 ig = 1, data%ng
            IF ( data%GXEQX( ig ) ) THEN
               f = f + data%GSCALE( ig ) * data%FT( ig )
            ELSE
               f = f + data%GSCALE( ig ) * data%GVALS( ig )
            END IF
  310    CONTINUE
      END IF
      IF ( GRAD ) THEN

!  Evaluate the group derivative values.

         IF ( .NOT. data%altriv ) CALL GROUP ( data%GVALS( 1 ), data%ng, &
               data%FT( 1 ), data%GPVALU( 1 ), icnt, &
               data%ITYPEG( 1 ), data%ISTGP( 1 ), &
               data%ICALCF( 1 ), &
               data%lcalcg, data%ng1, data%lcalcg, data%lcalcg, data%lgpvlu, &
               .TRUE., igstat )

!  Compute the gradient values. Initialize the gradient as zero.

         DO 120 j = 1, n
            G( j ) = zero
  120    CONTINUE   

!  Consider the IG-th group.

         DO 290 ig = 1, data%ng
            icon = data%KNDOFC( ig )

!  Consider only those groups in the objective function. 

            IF ( icon > 0 ) GO TO 290
            ig1 = ig + 1
            istrgv = data%IWORK( data%lstagv + ig )
            iendgv = data%IWORK( data%lstagv + ig1 ) - 1
            nelow = data%ISTADG( ig )
            nelup = data%ISTADG( ig1 ) - 1
            NONTRV = .NOT. data%GXEQX( ig )

!  Compute the first derivative of the group.

            gi = data%GSCALE( ig )
            IF ( NONTRV ) gi = gi  * data%GVALS( data%ng + ig ) 

!  The group has nonlinear elements.

            IF ( nelow <= nelup ) THEN
      CALL SETVI( iendgv - istrgv + 1, data%WRK( 1 ),  &
                            data%IWORK( data%lsvgrp + istrgv ), zero )

!  Loop over the group's nonlinear elements.

               DO 150 ii = nelow, nelup
                  iel = data%IELING( ii )
                  k = data%INTVAR( iel )
                  l = data%ISTAEV( iel )
                  nvarel = data%ISTAEV( iel + 1 ) - l
                  scalee = data%ESCALE( ii )
                  IF ( data%INTREP( iel ) ) THEN

!  The IEL-th element has an internal representation.

                     nin = data%INTVAR( iel + 1 ) - k
                     CALL RANGE ( iel, .TRUE., data%FUVALS( k ), &
                                  data%WRK( n + 1 ), nvarel, nin, &
                                  data%ITYPEE( iel ), &
                                  nin, nvarel )
!DIR$ IVDEP
                     DO 130 i = 1, nvarel
                        j = data%IELVAR( l )
                        data%WRK( j ) = data%WRK( j ) + &
                                           scalee * data%WRK( n + i )
                        l = l + 1
  130                CONTINUE
                  ELSE

!  The IEL-th element has no internal representation.

!DIR$ IVDEP
                     DO 140 i = 1, nvarel
                        j = data%IELVAR( l )
                        data%WRK( j ) = data%WRK( j ) + &
                                           scalee * data%FUVALS( k )
                        k = k + 1
                        l = l + 1
  140                CONTINUE
                  END IF
  150          CONTINUE

!  Include the contribution from the linear element.

!DIR$ IVDEP
               DO 160 k = data%ISTADA( ig ), &
                                  data%ISTADA( ig1 ) - 1
                  j = data%ICNA( k )
                  data%WRK( j ) = data%WRK( j ) + data%A( k )
  160          CONTINUE

!  Allocate a gradient.

!DIR$ IVDEP
               DO 190 i = istrgv, iendgv
                  ll = data%IWORK( data%lsvgrp + i )

!  Include the contributions from only the first n variables.

                  IF ( ll <= n ) &
                     G( ll ) = G( ll ) + gi * data%WRK( ll )
  190          CONTINUE

!  The group has only linear elements.

            ELSE

!  Allocate a gradient.

!DIR$ IVDEP
               DO 210 k = data%ISTADA( ig ), data%ISTADA( ig1 ) - 1
                  ll = data%ICNA( k )

!  Include the contributions from linear elements for only the first
!  n variables.

                  IF ( ll <= n ) &
                     G( ll ) = G( ll ) + gi * data%A( k )
  210          CONTINUE
            END IF
  290    CONTINUE
      ENDIF

!  Update the counters for the report tool.

      data%nc2of = data%nc2of + 1
      IF ( GRAD ) data%nc2og = data%nc2og + 1

      RETURN

!  Non-executable statements.

 2000 FORMAT( /  ' ** SUBROUTINE COFG: array length llogic too small.' &
              /  ' -- Minimization abandoned.' &
              /  ' -- Increase the parameter llogic by at least ', I8, &
                 ' and restart.' )

!  end of COFG.

      END

