! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE UOFG ( data, status, n, X, f, G, GRAD )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n
      INTEGER, INTENT( OUT ) :: status
      REAL ( KIND = wp ) :: f
      REAL ( KIND = wp ) :: X( n ), G( n )
      LOGICAL :: GRAD

!  Compute the value of the objective function and its gradient
!  for a function initially written in Standard Input Format (SIF).

!  G     is an array which gives the value of the gradient of
!        the objective function evaluated at X.
!        G(i) gives the partial derivative of the objective
!        function with respect to variable X(i).

!  Based on the subroutines ufn.f and ugr.f by Nick Gould, which are
!  in turn based on the subroutine SBMIN by Conn, Gould and Toint.

!  Ingrid Bongartz 
!  February 1993.

!  local variables.

      INTEGER :: i, j, ig, ifstat, igstat
      REAL ( KIND = wp ) :: ftt
      EXTERNAL :: RANGE 

!  There are non-trivial group functions.

      DO i = 1, MAX( data%nel, data%ng )
        data%ICALCF( i ) = i
      END DO

!  evaluate the element function values

      CALL ELFUN( data%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,          &
                  data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,          &
                  data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,           &
                  data%lelvar, data%lntvar, data%lstadh, data%lstep,           &
                  data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,          &
                  1, ifstat )

!  compute the group argument values ft

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
        data%FT( ig ) = ftt

!  record the derivatives of trivial groups

        IF ( data%GXEQX( ig ) ) data%GVALS( ig, 2 ) = 1.0_wp
      END DO

!  compute the group function values

!  all group functions are trivial

      IF ( data%altriv ) THEN
        f = DOT_PRODUCT( data%GSCALE( : data%ng ), data%FT( : data%ng ) )
        data%GVALS( : data%ng, 1 ) = data%FT( : data%ng )
        data%GVALS( : data%ng, 2 ) = 1.0_wp

!  evaluate the group function values

      ELSE
        CALL GROUP( data%GVALS, data%ng, data%FT, data%GPVALU, data%ng,        &
                    data%ITYPEG, data%ISTGP, data%ICALCF, data%ltypeg,         &
                    data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,         &
                    .FALSE., igstat )

!  compute the objective function value

        f = 0.0_wp
        DO ig = 1, data%ng
          IF ( data%GXEQX( ig ) ) THEN
            f = f + data%GSCALE( ig ) * data%FT( ig )
          ELSE
            f = f + data%GSCALE( ig ) * data%GVALS( ig, 1 )
          END IF
        END DO
      END IF

!  evaluate the element function derivatives

      IF ( grad ) THEN
        CALL ELFUN( data%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,        &
                    data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,        &
                    data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,         &
                    data%lelvar, data%lntvar, data%lstadh, data%lstep,         &
                    data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,        &
                    2, ifstat )

!  evaluate the group derivative values

        IF ( .NOT. data%altriv )                                               &
          CALL GROUP( data%GVALS, data%ng, data%FT, data%GPVALU, data%ng,      &
                      data%ITYPEG, data%ISTGP, data%ICALCF, data%ltypeg,       &
                      data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,       &
                      .TRUE., igstat )

!  compute the gradient values

        CALL ELGRD( n, data%ng, data%firstg, data%ICNA( 1 ), data%licna, &
                      data%ISTADA( 1 ), data%lstada, data%IELING( 1 ), &
                      data%leling, data%ISTADG( 1 ), data%lstadg, &
                      data%ITYPEE( 1 ), data%lintre, &
                      data%ISTAEV( 1 ), data%lstaev, &
                      data%IELVAR( 1 ), data%lelvar, &
                      data%INTVAR( 1 ), data%lntvar, &
                      data%IWORK(LSVGRP + 1), data%lnvgrp, &
                      data%IWORK( data%lstajc + 1 ), data%lnstjc, &
                      data%IWORK( data%lstagv + 1 ), data%lnstgv, &
                      data%A( 1 ), data%la, &
                      data%GVALS( : , 2 ), data%lgvals, &
                      data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ), &
                      data%GSCALE( 1 ), data%lgscal, &
                      data%ESCALE( 1 ), data%lescal, &
                      data%FUVALS( data%lgrjac + 1 ), &
                      data%lngrjc, data%WRK( 1 ), &
                      data%WRK( n + 1 ), data%maxsel, &
                      data%GXEQX( 1 ), data%lgxeqx, &
                      data%INTREP( 1 ), data%lintre, RANGE )
        data%firstg = .FALSE.

!  store the gradient value

        DO i = 1, n
          G( i ) = data%FUVALS( data%lggfx + i )
        END DO
      END IF

!  update the counters for the report tool

      data%nc2of = data%nc2of + 1
      IF ( grad ) data%nc2og = data%nc2og + 1
      status = 0
      RETURN

!  end of UOFG

      END
