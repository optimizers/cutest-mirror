! ( Last modified on 10 Sepc 2004 at 16:55:38 )
!  Correction: 10/Sep/2004: undeclared integers variables declared
      SUBROUTINE CJPROD( data, n, m, GOTJ, JTRANS, X, V, lv, R, lr )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n, m, lv, lr
      LOGICAL :: GOTJ, JTRANS
      REAL ( KIND = wp ) :: X( n ), V( lv ), R( lr )

!  Compute the matrix-vector product between the Jacobian matrix
!  of the constraints (JTRANS = .FALSE.), or its transpose 
! (JTRANS = .TRUE.) for the problem, and a given vector P. 
!  The result is placed in R. If GOTJ is .TRUE. the first derivatives 
!  are assumed to have already been computed. If the user is unsure, 
!  set GOTJ = .FALSE. the first time a product is required with the 
!  Jacobian evaluated at X. X is not used if GOTJ = .TRUE.

!  Nick Gould, for GOT/CUTEr productions.
!  June, 2003.


! ---------------------------------------------------------------------




! ---------------------------------------------------------------------



! ---------------------------------------------------------------------


! ---------------------------------------------------------------------





!  Integer variables from the PRFCTS common block.


!  Local variables

      INTEGER :: i, ig, j, icon, k, ig1, ii
      INTEGER :: l, iel, nvarel, nin
      INTEGER :: ifstat, igstat
      REAL ( KIND = wp ) :: zero, one, ftt, prod, scalee
      PARAMETER ( zero = 0.0_wp, one = 1.0_wp )
      EXTERNAL :: RANGE 
      IF ( data%numcon == 0 ) RETURN

!  Check input data.

      IF ( ( JTRANS .AND. lv < m ) .OR.  &
 ( .NOT. JTRANS .AND. lv < n ) ) THEN
         IF ( iout > 0 ) WRITE( iout, 2010 )
         STOP
      END IF
      IF ( ( JTRANS .AND. lr < n ) .OR.  &
 ( .NOT. JTRANS .AND. lr < m ) ) THEN
         IF ( iout > 0 ) WRITE( iout, 2020 )
         STOP
      END IF

!  There are non-trivial group functions.

      IF ( .NOT. GOTJ ) THEN
         DO 10 i = 1, MAX( data%nel, data%ng )
           data%ICALCF( i ) = i
   10    CONTINUE

!  Evaluate the element function values.

         CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), data%nel, &
                      data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                      data%IELVAR( 1 ), data%INTVAR( 1 ), &
                      data%ISTADH( 1 ), data%ISTEP( 1 ), &
                      data%ICALCF( 1 ),  &
                      data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                      data%lntvar, data%lintre, lfuval, data%lvscal, data%lepvlu,  &
                      1, ifstat )

!  Evaluate the element function values.

         CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), data%nel, &
                      data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                      data%IELVAR( 1 ), data%INTVAR( 1 ), &
                      data%ISTADH( 1 ), data%ISTEP( 1 ), &
                      data%ICALCF( 1 ),  &
                      data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                      data%lntvar, data%lintre, lfuval, data%lvscal, data%lepvlu,  &
                      3, ifstat )

!  Compute the group argument values ft.

         DO 70 ig = 1, data%ng
            ftt = - data%B( ig )

!  Include the contribution from the linear element.

            DO 30 j = data%ISTADA( ig ), data%ISTADA( ig + 1 ) - 1
               ftt = ftt + data%A( j ) * X( data%ICNA( j ) )
   30       CONTINUE

!  Include the contributions from the nonlinear elements.

            DO 60 j = data%ISTADG( ig ), data%ISTADG( ig + 1 ) - 1
               ftt = ftt + data%ESCALE( j ) * data%FUVALS( data%IELING( J))
   60       CONTINUE
            data%FT( ig ) = ftt

!  Record the derivatives of trivial groups.

            IF ( data%GXEQX( ig ) ) THEN
               data%GVALS( data%ng + ig ) = one
            END IF
   70    CONTINUE

!  Evaluate the group derivative values.

         IF ( .NOT. data%altriv ) CALL GROUP ( data%GVALS( 1 ), data%ng, &
               data%FT( 1 ), data%GPVALU( 1 ), data%ng, &
               data%ITYPEG( 1 ), data%ISTGP( 1 ), &
               data%ICALCF( 1 ), &
               data%lcalcg, data%ng1, data%lcalcg, data%lcalcg, data%lgpvlu, &
               .TRUE., igstat )
      END IF

!  Ensure that there is sufficient space.

      IF ( data%lwk2 < n ) THEN
         IF ( iout > 0 ) WRITE( iout, 2000 )
         STOP
      END IF

!  Form the product r = J(transpose) v

      IF ( JTRANS ) THEN

!  Initialize R

         DO 110 i = 1, n
            R( i ) = zero
  110    CONTINUE

!  Consider the IG-th group.

         DO 190 ig = 1, data%ng
            icon = data%KNDOFC( ig )
            IF ( icon > 0 ) THEN
               ig1 = ig + 1


!  Compute the product of v(i) with the (scaled) group derivative

               IF ( data%GXEQX( ig ) ) THEN
                  prod = V( icon ) * data%GSCALE( ig )
               ELSE
                  prod = V( icon ) * data%GSCALE( ig ) * data%GVALS( data%ng + ig )
               END IF

!  Loop over the group's nonlinear elements.

               DO 150 ii = data%ISTADG( ig ),  &
                           data%ISTADG( ig1 ) - 1
                  iel = data%IELING( ii )
                  k = data%INTVAR( iel )
                  l = data%ISTAEV( iel )
                  nvarel = data%ISTAEV( iel + 1 ) - l
                  scalee = data%ESCALE( ii ) * prod
                  IF ( data%INTREP( iel ) ) THEN

!  The IEL-th element has an internal representation.

                     nin = data%INTVAR( iel + 1 ) - k
                     CALL RANGE ( iel, .TRUE., data%FUVALS( k ), &
                                  data%WRK( 1 ), nvarel, nin, &
                                  data%ITYPEE( iel ), &
                                  nin, nvarel )
!DIR$ IVDEP
                     DO 130 i = 1, nvarel
                        j = data%IELVAR( l )
                        R( j ) = R( j ) + scalee * data%WRK( i )
                        l = l + 1
  130                CONTINUE
                  ELSE

!  The IEL-th element has no internal representation.

!DIR$ IVDEP
                     DO 140 i = 1, nvarel
                        j = data%IELVAR( l )
                        R( j ) = R( j ) + scalee * data%FUVALS( k )
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
                  R( j ) = R( j ) + data%A( k ) * prod
  160          CONTINUE
            END IF
  190    CONTINUE

!  Form the product r = j v

      ELSE

!  Consider the IG-th group.

         DO 290 ig = 1, data%ng
            icon = data%KNDOFC( ig )
            IF ( icon > 0 ) THEN
               ig1 = ig + 1
               prod = zero

!  Compute the first derivative of the group.

!  Loop over the group's nonlinear elements.

               DO 250 ii = data%ISTADG( ig ),  &
                           data%ISTADG( ig1 ) - 1
                  iel = data%IELING( ii )
                  k = data%INTVAR( iel )
                  l = data%ISTAEV( iel )
                  nvarel = data%ISTAEV( iel + 1 ) - l
                  scalee = data%ESCALE( ii )
                  IF ( data%INTREP( iel ) ) THEN

!  The IEL-th element has an internal representation.

                     nin = data%INTVAR( iel + 1 ) - k
                     CALL RANGE ( iel, .TRUE., data%FUVALS( k ), &
                                  data%WRK( 1 ), nvarel, nin, &
                                  data%ITYPEE( iel ), &
                                  nin, nvarel )
!DIR$ IVDEP
                     DO 230 i = 1, nvarel
                        prod = prod + V( data%IELVAR( l ) ) * &
                                       scalee * data%WRK( i )
                        l = l + 1
  230                CONTINUE
                  ELSE

!  The IEL-th element has no internal representation.

!DIR$ IVDEP
                     DO 240 i = 1, nvarel
                        prod = prod + V( data%IELVAR( l ) ) * &
                                   scalee * data%FUVALS( k )
                        k = k + 1
                        l = l + 1
  240                CONTINUE
                  END IF
  250          CONTINUE

!  Include the contribution from the linear element.

!DIR$ IVDEP
               DO 260 k = data%ISTADA( ig ), &
                          data%ISTADA( ig1 ) - 1
                  prod = prod + V( data%ICNA( k ) ) * data%A( k )
  260          CONTINUE

!  Multiply the product by the (scaled) group derivative

               IF ( data%GXEQX( ig ) ) THEN
                  R( icon ) = prod * data%GSCALE( ig )
               ELSE
                  R( icon ) = prod * data%GSCALE( ig ) * data%GVALS( data%ng + ig )
               END IF
            END IF
  290    CONTINUE
      END IF

!  Update the counters for the report tool.

      IF ( .NOT. GOTJ ) THEN
         data%nc2og = data%nc2og + 1
         data%nc2cg = data%nc2cg + data%pnc
      END IF
      RETURN

! Non-executable statements.

 2000 FORMAT( ' ** SUBROUTINE CJPROD: Increase the size of WK ' )
 2010 FORMAT( ' ** SUBROUTINE CJPROD: Increase the size of V ' )
 2020 FORMAT( ' ** SUBROUTINE CJPROD: Increase the size of R ' )

!  end of CJPROD.

      END
