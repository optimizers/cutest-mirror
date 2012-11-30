! THIS VERSION: CUTEST 1.0 - 28/11/2012 AT 08:00 GMT.

!-*-*-*-*-*-*-*-  C U T E S T    C O F G    S U B R O U T I N E  -*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal authors: Ingrid Bongartz and Nick Gould

!  History -
!   fortran 77 version originally released in CUTE, April 1992
!   fortran 2003 version released in CUTEst, 28th November 2012

      SUBROUTINE COFG( data, status, n, X, f, G, grad )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      TYPE ( CUTEST_data_type ), INTENT( INOUT ) :: data
      INTEGER, INTENT( IN ) :: n
      INTEGER, INTENT( OUT ) :: status
      REAL ( KIND = wp ), INTENT( OUT ) :: f
      LOGICAL, INTENT( IN ) :: grad
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: X
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: G

!  ---------------------------------------------------------------
!  compute the value of the objective function and its gradient
!  for a function initially written in Standard Input Format (SIF)

!  G is an array which gives the value of the gradient of the
!    objective function evaluated at X. G(i) gives the partial 
!    derivative of the objective function wrt variable X(i)
!  ---------------------------------------------------------------

!  local variables

      INTEGER :: i, j, iel, k, ig, ii, ig1, l, ll, icon, icnt, ifstat, igstat
      INTEGER :: nin, nvarel, nelow, nelup, istrgv, iendgv
      EXTERNAL :: RANGE 
      REAL ( KIND = wp ) :: ftt, gi, scalee

!  identify which elements are included in objective function. Use LOGIC
!  to keep track of elements already included

      data%LOGIC( : data%nel ) = .FALSE.

!  now identify elements in objective function groups

      icnt = 0
      DO ig = 1, data%ng
        IF ( data%KNDOFC( ig ) == 0 ) THEN
          nelow = data%ISTADG( ig )
          nelup = data%ISTADG( ig + 1 ) - 1
          DO ii = nelow, nelup
            iel = data%IELING( ii )
            IF ( .NOT. data%LOGIC( iel ) ) THEN
              data%LOGIC( iel ) = .TRUE.
              icnt = icnt + 1
              data%ICALCF( icnt ) = iel
            END IF
          END DO
        END IF
      END DO

!  evaluate the element function values

      CALL ELFUN( data%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,          &
                  data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,          &
                  data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,           &
                  data%lelvar, data%lntvar, data%lstadh, data%lstep,           &
                  data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,          &
                  1, ifstat )
      IF ( ifstat /= 0 ) GO TO 930

!  evaluate the element function derivatives

      IF ( GRAD )                                                              &
        CALL ELFUN( data%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,        &
                    data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,        &
                    data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,         &
                    data%lelvar, data%lntvar, data%lstadh, data%lstep,         &
                    data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,        &
                    2, ifstat )
      IF ( ifstat /= 0 ) GO TO 930

!  compute the group argument values ft

      IF ( data%numcon > 0 ) THEN
        DO ig = 1, data%ng
          ftt = 0.0_wp

!  consider only those groups in the objective function

          IF ( data%KNDOFC( ig ) == 0 ) THEN
            ftt = - data%B( ig )

!  include the contribution from the linear element only if the variable 
!  belongs to the first n variables

            DO i = data%ISTADA( ig ), data%ISTADA( ig + 1 ) - 1
              j = data%ICNA( i ) 
              IF ( j <= n ) ftt = ftt + data%A( i ) * X( j )
            END DO

!  include the contributions from the nonlinear elements

            DO i = data%ISTADG( ig ), data%ISTADG( ig + 1 ) - 1
              ftt = ftt + data%ESCALE( i ) * data%FUVALS( data%IELING( i ) )
            END DO

!  record the derivatives of trivial groups

            IF ( data%GXEQX( ig ) ) data%GVALS( ig, 2 ) = 1.0_wp
          END IF
          data%FT( ig ) = ftt
        END DO
      ELSE

!  there are no constraints, so we need not check data%KNDOFC(ig)

        DO ig = 1, data%ng
          ftt = - data%B( ig )

!  include the contribution from the linear element only if the variable 
!  belongs to the first n variables

          DO i = data%ISTADA( ig ), data%ISTADA( ig + 1 ) - 1
            j = data%ICNA( i )
            IF ( j <= n ) ftt = ftt + data%A( i ) * X( j )
          END DO

!  include the contributions from the nonlinear elements

          DO i = data%ISTADG( ig ), data%ISTADG( ig + 1 ) - 1
            ftt = ftt + data%ESCALE( i ) * data%FUVALS( data%IELING( i ) )
          END DO

!  record the derivatives of trivial groups

          IF ( data%GXEQX( ig ) ) data%GVALS( ig, 2 ) = 1.0_wp
          data%FT( ig ) = ftt
        END DO
      END IF

!  compute the group function values

!  all group functions are trivial

      IF ( data%altriv ) THEN
        data%GVALS( : data%ng, 1 ) = data%FT( : data%ng )
        data%GVALS( : data%ng, 2 ) = 1.0_wp
      ELSE

!  evaluate the group function values. Only evaluate groups belonging to the 
!  objective function

        icnt = 0
        DO ig = 1, data%ng
          IF ( data%KNDOFC( ig ) == 0 ) THEN
            icnt = icnt + 1
            data%ICALCF( icnt ) = ig
          END IF 
        END DO
        CALL GROUP( data%GVALS, data%ng, data%FT, data%GPVALU, data%ng,        &
                    data%ITYPEG, data%ISTGP, data%ICALCF, data%ltypeg,         &
                    data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,         &
                    .FALSE., igstat )
        IF ( igstat /= 0 ) GO TO 930
      END IF

!  compute the objective function value

      f = 0.0_wp

!  there are constraints

      IF ( data%numcon > 0 ) THEN
        DO ig = 1, data%ng
          IF ( data%KNDOFC( ig ) == 0 ) THEN
            IF ( data%GXEQX( ig ) ) THEN
              f = f + data%GSCALE( ig ) * data%FT( ig )
            ELSE
              f = f + data%GSCALE( ig ) * data%GVALS( ig, 1 )
            END IF
          END IF 
        END DO
      ELSE

!  there are no constraints, so we need not check data%KNDOFC( ig )

        DO ig = 1, data%ng
          IF ( data%GXEQX( ig ) ) THEN
            f = f + data%GSCALE( ig ) * data%FT( ig )
          ELSE
            f = f + data%GSCALE( ig ) * data%GVALS( ig, 1 )
          END IF
        END DO
      END IF

!  evaluate the group derivative values

      IF ( grad ) THEN
        IF ( .NOT. data%altriv ) THEN
          CALL GROUP( data%GVALS, data%ng, data%FT, data%GPVALU, data%ng,      &
                      data%ITYPEG, data%ISTGP, data%ICALCF, data%ltypeg,       &
                      data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,       &
                      .TRUE., igstat )
          IF ( igstat /= 0 ) GO TO 930
        END IF

!  compute the gradient values. Initialize the gradient as zero

        G( : n ) = 0.0_wp

!  consider the IG-th group

        DO ig = 1, data%ng
          icon = data%KNDOFC( ig )

!  consider only those groups in the objective function

          IF ( icon > 0 ) CYCLE
          ig1 = ig + 1
          istrgv = data%ISTAGV( ig )
          iendgv = data%ISTAGV( ig1 ) - 1
          nelow = data%ISTADG( ig )
          nelup = data%ISTADG( ig1 ) - 1

!  compute the first derivative of the group

          gi = data%GSCALE( ig )
          IF ( .NOT. data%GXEQX( ig ) ) gi = gi  * data%GVALS( ig, 2 ) 

!  the group has nonlinear elements

          IF ( nelow <= nelup ) THEN
            data%W_ws( data%ISVGRP( istrgv : iendgv ) ) = 0.0_wp

!  loop over the group's nonlinear elements

            DO ii = nelow, nelup
              iel = data%IELING( ii )
              k = data%INTVAR( iel )
              l = data%ISTAEV( iel )
              nvarel = data%ISTAEV( iel + 1 ) - l
              scalee = data%ESCALE( ii )
              IF ( data%INTREP( iel ) ) THEN

!  the iel-th element has an internal representation

                nin = data%INTVAR( iel + 1 ) - k
                CALL RANGE( iel, .TRUE., data%FUVALS( k ), data%W_el,          &
                            nvarel, nin, data%ITYPEE( iel ), nin, nvarel )
!DIR$ IVDEP
                DO i = 1, nvarel
                  j = data%IELVAR( l )
                  data%W_ws( j ) = data%W_ws( j ) + scalee * data%W_el( i )
                  l = l + 1
                END DO
              ELSE

!  the iel-th element has no internal representation

!DIR$ IVDEP
                DO i = 1, nvarel
                  j = data%IELVAR( l )
                  data%W_ws( j ) = data%W_ws( j ) + scalee * data%FUVALS( k )
                   k = k + 1 ; l = l + 1
                END DO
              END IF
            END DO

!  include the contribution from the linear element

!DIR$ IVDEP
            DO k = data%ISTADA( ig ), data%ISTADA( ig1 ) - 1
             j = data%ICNA( k )
             data%W_ws( j ) = data%W_ws( j ) + data%A( k )
            END DO

!  allocate a gradient

!DIR$ IVDEP
            DO i = istrgv, iendgv
              ll = data%ISVGRP( i )

!  include the contributions from only the first n variables

              IF ( ll <= n ) G( ll ) = G( ll ) + gi * data%W_ws( ll )
            END DO

!  the group has only linear elements

          ELSE

!  allocate a gradient

!DIR$ IVDEP
            DO k = data%ISTADA( ig ), data%ISTADA( ig1 ) - 1
              ll = data%ICNA( k )

!  include the contributions from linear elements for only the first n 
!  variables

              IF ( ll <= n ) G( ll ) = G( ll ) + gi * data%A( k )
            END DO
          END IF
        END DO
      END IF

!  update the counters for the report tool

      data%nc2of = data%nc2of + 1
      IF ( grad ) data%nc2og = data%nc2og + 1

      status = 0
      RETURN

!  unsuccessful returns

  930 CONTINUE
      IF ( data%out > 0 ) WRITE( data%out,                                     &
        "( ' ** SUBROUTINE COFG: error flag raised during SIF evaluation' )" )
      status = 3
      RETURN

!  end of subroutine COFG

      END SUBROUTINE COFG

