! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CCIFG( data, status, n, icon, X, ci, GCI, grad )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n, icon
      INTEGER, INTENT( OUT ) :: status
      REAL ( KIND = wp ) :: ci
      REAL ( KIND = wp ) :: X( n ), GCI( n )
      LOGICAL :: grad

!  Evaluate constraint function icon and possibly its gradient,
!  for constraints initially written in Standard Input Format (SIF).  
!  The constraint gradient is stored as a dense vector in array GCI;
!  that is, GCI( j ) is the partial derivative of constraint icon
!  with respect to variable j.
! (Subroutine CSCIFG performs the same calculations for a sparse
!  constraint gradient vector.)

!  Ingrid Bongartz
!  September 1994.

!  local variables.

      INTEGER :: i, j, iel, k, ig, ii, ig1, l, ll, neling
      INTEGER :: nin, nvarel, nelow, nelup, istrgv, iendgv
      INTEGER :: ifstat, igstat
      LOGICAL :: NONTRV
      EXTERNAL :: RANGE 
      REAL ( KIND = wp ) :: ftt, gi, scalee

!  Return if there are no constraints.

      IF ( data%numcon == 0 ) RETURN

!  check input parameters

      IF ( icon <= 0 ) THEN
        IF ( data%out > 0 ) WRITE( data%out, 2000 )
        status = 2 ; RETURN
      END IF

!  find group index ig of constraint icon

      ig = 0
      DO i = 1, data%ng
        IF ( data%KNDOFC( i ) == icon ) THEN
          ig = i
          EXIT
        END IF
      END DO
      IF ( ig == 0 ) THEN
        IF ( data%out > 0 ) WRITE( data%out, 2000 )
        status = 2 ; RETURN
      END IF

!  determine nonlinear elements in group ig. Record their indices in 
!  data%IWORK( ICALCF )

      nelow = data%ISTADG( ig )
      nelup = data%ISTADG( ig + 1 ) - 1
      neling = nelup - nelow + 1
      j = nelow - 1
      DO i = 1, neling
        j = j + 1
        data%ICALCF( i ) = data%IELING( j )
      END DO

!  evaluate the element functions

      CALL ELFUN( data%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,          &
                  data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,          &
                  data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,           &
                  data%lelvar, data%lntvar, data%lstadh, data%lstep,           &
                  data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,          &
                  1, ifstat )
      IF ( ifstat /= 0 ) GO TO 930

!  compute the group argument value FTT. Consider only the group associated 
!  with constraint icon

      ftt = - data%B( ig )

!  include the contribution from the linear element only if the variable 
!  belongs to the first n variables

      DO i = data%ISTADA( ig ), data%ISTADA( ig + 1 ) - 1
        j = data%ICNA( i )
        IF ( j <= n ) ftt = ftt + data%A( i ) * X( j )
      END DO

!  Include the contributions from the nonlinear elements.

      DO i = nelow, nelup
         ftt = ftt + data%ESCALE( i ) * data%FUVALS( data%IELING( i ) )
      END DO
      data%FT( ig ) = ftt

!  If ig is a trivial group, record the function value and derivative.

      IF ( data%GXEQX( ig ) ) THEN
        data%GVALS( ig, 1 ) = data%FT( ig )
        data%GVALS( ig, 2 ) = 1.0_wp

!  otherwise, evaluate group ig

      ELSE
        CALL GROUP( data%GVALS, data%ng, data%FT, data%GPVALU, data%ng,        &
                    data%ITYPEG, data%ISTGP, data%ICALCF, data%ltypeg,         &
                    data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,         &
                    .FALSE., igstat )
        IF ( igstat /= 0 ) GO TO 930
      END IF

!  Compute the constraint function value.

      IF ( data%GXEQX( ig ) ) THEN
        ci = data%GSCALE( ig ) * data%FT( ig )
      ELSE
        ci = data%GSCALE( ig ) * data%GVALS( ig, 1 )

!  Update the constraint function evaluation counter

        data%nc2cf = data%nc2cf + 1
      END IF
      IF ( grad ) THEN

!  Update the constraint gradient evaluation counter

        data%nc2cg = data%nc2cg + 1

!  evaluate the element function derivatives

        CALL ELFUN( data%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,        &
                    data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,        &
                    data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,         &
                    data%lelvar, data%lntvar, data%lstadh, data%lstep,         &
                    data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,        &
                    2, ifstat )
        IF ( ifstat /= 0 ) GO TO 930

!  evaluate the group derivative

        IF ( .NOT. data%GXEQX( ig ) ) THEN
          CALL GROUP( data%GVALS, data%ng, data%FT, data%GPVALU, data%ng,      &
                      data%ITYPEG, data%ISTGP, data%ICALCF, data%ltypeg,       &
                      data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,       &
                      .TRUE., igstat )
          IF ( igstat /= 0 ) GO TO 930
        END IF

!  compute the gradient. Initialize the gradient vector as zero

        GCI( : n ) = 0.0_wp

!  consider only group ig

        ig1 = ig + 1
        istrgv = data%IWORK( data%lstagv + ig )
        iendgv = data%IWORK( data%lstagv + ig1 ) - 1
        nontrv = .NOT. data%GXEQX( ig )

!  compute the first derivative of the group

        gi = data%GSCALE( ig )
        IF ( nontrv ) gi = gi  * data%GVALS( ig, 2 )

!  the group has nonlinear elements

        IF ( nelow <= nelup ) THEN
          data%WRK( data%IWORK( data%lsvgrp + istrgv :                        &
                                data%lsvgrp + iendgv ) ) = 0.0_wp

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
              CALL RANGE( iel, .TRUE., data%FUVALS( k ),                       &
                           data%WRK( n + 1 ), nvarel, nin,                     &
                           data%ITYPEE( iel ), nin, nvarel )
!DIR$ IVDEP
              DO i = 1, nvarel
                j = data%IELVAR( l )
                data%WRK( j ) = data%WRK( j ) + scalee * data%WRK( n + i )
                l = l + 1
              END DO

!  the iel-th element has no internal representation

            ELSE
!DIR$ IVDEP
              DO i = 1, nvarel
                j = data%IELVAR( l )
                data%WRK( j ) = data%WRK( j ) + scalee * data%FUVALS( k )
                k = k + 1 ; l = l + 1
              END DO
            END IF
          END DO

!  include the contribution from the linear element

!DIR$ IVDEP
          DO k = data%ISTADA( ig ), data%ISTADA( ig1 ) - 1
            j = data%ICNA( k )
            data%WRK( j ) = data%WRK( j ) + data%A( k )
          END DO

!  allocate a gradient

!DIR$ IVDEP
          DO i = istrgv, iendgv
            ll = data%IWORK( data%lsvgrp + i )

!  include contributions from the first n variables only

            IF ( ll <= n ) GCI( ll ) = gi * data%WRK( ll )
          END DO

!  the group has only linear elements

        ELSE

!  allocate a gradient

!DIR$ IVDEP
          DO k = data%ISTADA( ig ), data%ISTADA( ig1 ) - 1
            ll = data%ICNA( k )

!  include contributions from the first n variables only

            IF ( ll <= n ) GCI( ll ) = gi * data%A( k )
          END DO
        END IF
      END IF
      status = 0
      RETURN

!  unsuccessful returns

  930 CONTINUE
      IF ( data%out > 0 ) WRITE( data%out,                                     &
        "( ' ** SUBROUTINE CCIFG: error flag raised during SIF evaluation' )" )
      status = 3
      RETURN

!  non-executable statements

 2000 FORMAT( ' ** SUBROUTINE CCIFG: invalid constraint index icon ' )

!  end of subroutine CCIFG

      END SUBROUTINE CCIFG
