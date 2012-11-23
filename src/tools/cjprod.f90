! ( Last modified on 10 Sepc 2004 at 16:55:38 )
!  Correction: 10/Sep/2004: undeclared integers variables declared
      SUBROUTINE CJPROD( data, status, n, m, gotj, jtrans, X, V, lv, R, lr )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n, m, lv, lr
      INTEGER, INTENT( OUT ) :: status
      LOGICAL :: gotj, jtrans
      REAL ( KIND = wp ) :: X( n ), V( lv ), R( lr )

!  Compute the matrix-vector product between the Jacobian matrix
!  of the constraints (JTRANS = .FALSE.), or its transpose 
! (JTRANS = .TRUE.) for the problem, and a given vector P. 
!  The result is placed in R. If gotj is .TRUE. the first derivatives 
!  are assumed to have already been computed. If the user is unsure, 
!  set gotj = .FALSE. the first time a product is required with the 
!  Jacobian evaluated at X. X is not used if gotj = .TRUE.

!  Nick Gould, for GOT/CUTEr productions.
!  June, 2003.

!  Integer variables from the PRFCTS common block.


!  Local variables

      INTEGER :: i, ig, j, icon, k, ig1, ii
      INTEGER :: l, iel, nvarel, nin
      INTEGER :: ifstat, igstat
      REAL ( KIND = wp ) :: ftt, prod, scalee
      EXTERNAL :: RANGE 
      IF ( data%numcon == 0 ) RETURN

!  check input data

      IF ( ( jtrans .AND. lv < m ) .OR. ( .NOT. jtrans .AND. lv < n ) ) THEN
         IF ( data%out > 0 ) WRITE( data%out, 2010 )
         status = 2 ; RETURN
      END IF
      IF ( ( jtrans .AND. lr < n ) .OR.( .NOT. jtrans .AND. lr < m ) ) THEN
         IF ( data%out > 0 ) WRITE( data%out, 2020 )
         status = 2 ; RETURN
      END IF

!  there are non-trivial group functions

      IF ( .NOT. gotj ) THEN
         DO i = 1, MAX( data%nel, data%ng )
           data%ICALCF( i ) = i
         END DO

!  evaluate the element function values

        CALL ELFUN( data%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,        &
                    data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,        &
                    data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,         &
                    data%lelvar, data%lntvar, data%lstadh, data%lstep,         &
                    data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,        &
                    1, ifstat )
      IF ( ifstat /= 0 ) GO TO 930

!  evaluate the element function values

        CALL ELFUN( data%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,        &
                    data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,        &
                    data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,         &
                    data%lelvar, data%lntvar, data%lstadh, data%lstep,         &
                    data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,        &
                    3, ifstat )
      IF ( ifstat /= 0 ) GO TO 930

!  compute the group argument values ft

        DO ig = 1, data%ng
          ftt = - data%B( ig )

!  include the contribution from the linear element

          DO j = data%ISTADA( ig ), data%ISTADA( ig + 1 ) - 1
            ftt = ftt + data%A( j ) * X( data%ICNA( j ) )
          END DO

!  include the contributions from the nonlinear elements

          DO j = data%ISTADG( ig ), data%ISTADG( ig + 1 ) - 1
            ftt = ftt + data%ESCALE( j ) * data%FUVALS( data%IELING( J ) )
          END DO
          data%FT( ig ) = ftt

!  record the derivatives of trivial groups

          IF ( data%GXEQX( ig ) ) data%GVALS( ig, 2 ) = 1.0_wp
        END DO

!  evaluate the group derivative values

        IF ( .NOT. data%altriv ) THEN
          CALL GROUP( data%GVALS, data%ng, data%FT, data%GPVALU, data%ng,      &
                      data%ITYPEG, data%ISTGP, data%ICALCF, data%ltypeg,       &
                      data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,       &
                      .TRUE., igstat )
         IF ( igstat /= 0 ) GO TO 930
        END IF
      END IF

!  ensure that there is sufficient space

      IF ( data%lwk2 < n ) THEN
        IF ( data%out > 0 ) WRITE( data%out, 2000 )
        status = 2 ; RETURN
      END IF

!  form the product r = J(transpose) v

      IF ( jtrans ) THEN

!  initialize R

        R( : n ) = 0.0_wp

!  consider the ig-th group

        DO ig = 1, data%ng
          icon = data%KNDOFC( ig )
          IF ( icon > 0 ) THEN
            ig1 = ig + 1

!  compute the product of v(i) with the (scaled) group derivative

            IF ( data%GXEQX( ig ) ) THEN
              prod = V( icon ) * data%GSCALE( ig )
            ELSE
              prod = V( icon ) * data%GSCALE( ig ) * data%GVALS( ig, 2 )
            END IF

!  loop over the group's nonlinear elements

            DO ii = data%ISTADG( ig ), data%ISTADG( ig1 ) - 1
              iel = data%IELING( ii )
              k = data%INTVAR( iel )
              l = data%ISTAEV( iel )
              nvarel = data%ISTAEV( iel + 1 ) - l
              scalee = data%ESCALE( ii ) * prod
              IF ( data%INTREP( iel ) ) THEN

!  the iel-th element has an internal representation

                nin = data%INTVAR( iel + 1 ) - k
                CALL RANGE ( iel, .TRUE., data%FUVALS( k ),                    &
                             data%WRK( 1 ), nvarel, nin,                       &
                             data%ITYPEE( iel ), nin, nvarel )
!DIR$ IVDEP
                DO i = 1, nvarel
                  j = data%IELVAR( l )
                  R( j ) = R( j ) + scalee * data%WRK( i )
                  l = l + 1
                END DO
              ELSE

!  the iel-th element has no internal representation

!DIR$ IVDEP
                DO i = 1, nvarel
                  j = data%IELVAR( l )
                  R( j ) = R( j ) + scalee * data%FUVALS( k )
                  k = k + 1 ; l = l + 1
                 END DO
              END IF
            END DO

!  include the contribution from the linear element

!DIR$ IVDEP
            DO k = data%ISTADA( ig ), data%ISTADA( ig1 ) - 1
              j = data%ICNA( k )
              R( j ) = R( j ) + data%A( k ) * prod
            END DO
          END IF
        END DO

!  Form the product r = j v

      ELSE

!  consider the IG-th group

        DO ig = 1, data%ng
          icon = data%KNDOFC( ig )
          IF ( icon > 0 ) THEN
            ig1 = ig + 1
            prod = 0.0_wp

!  compute the first derivative of the group

!  loop over the group's nonlinear elements

            DO ii = data%ISTADG( ig ), data%ISTADG( ig1 ) - 1
              iel = data%IELING( ii )
              k = data%INTVAR( iel )
              l = data%ISTAEV( iel )
              nvarel = data%ISTAEV( iel + 1 ) - l
              scalee = data%ESCALE( ii )
              IF ( data%INTREP( iel ) ) THEN

!  the iel-th element has an internal representation

                nin = data%INTVAR( iel + 1 ) - k
                CALL RANGE ( iel, .TRUE., data%FUVALS( k ),                 &
                             data%WRK( 1 ), nvarel, nin,                    &
                             data%ITYPEE( iel ), nin, nvarel )
!DIR$ IVDEP
                DO i = 1, nvarel
                  prod = prod + V( data%IELVAR( l ) ) * scalee * data%WRK( i )
                  l = l + 1
                END DO
              ELSE

!  the iel-th element has no internal representation

!DIR$ IVDEP
                DO i = 1, nvarel
                  prod = prod + V( data%IELVAR( l ) ) * scalee * data%FUVALS( k)
                  k = k + 1
                  l = l + 1
                END DO
              END IF
            END DO

!  include the contribution from the linear element

!DIR$ IVDEP
            DO k = data%ISTADA( ig ), data%ISTADA( ig1 ) - 1
              prod = prod + V( data%ICNA( k ) ) * data%A( k )
            END DO

!  multiply the product by the (scaled) group derivative

            IF ( data%GXEQX( ig ) ) THEN
               R( icon ) = prod * data%GSCALE( ig )
            ELSE
               R( icon ) = prod * data%GSCALE( ig ) * data%GVALS( ig, 2 )
            END IF
          END IF
        END DO
      END IF

!  update the counters for the report tool

      IF ( .NOT. gotj ) THEN
         data%nc2og = data%nc2og + 1
         data%nc2cg = data%nc2cg + data%pnc
      END IF
      status = 0
      RETURN

!  unsuccessful returns

  930 CONTINUE
      IF ( data%out > 0 ) WRITE( data%out,                                     &
        "( ' ** SUBROUTINE CJPROD: error flag raised during SIF evaluation' )" )
      status = 3
      RETURN

! Non-executable statements.

 2000 FORMAT( ' ** SUBROUTINE CJPROD: Increase the size of WK' )
 2010 FORMAT( ' ** SUBROUTINE CJPROD: Increase the size of V' )
 2020 FORMAT( ' ** SUBROUTINE CJPROD: Increase the size of R' )

!  end of subroutine CJPROD

      END SUBROUTINE CJPROD
