! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CCFSG( data, status, n, m, X, C, nnzj, lcjac, CJAC,           &
                        INDVAR, INDFUN, grad )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n, m, nnzj, lcjac
      INTEGER, INTENT( OUT ) :: status
      LOGICAL :: grad
      INTEGER :: INDVAR( lcjac ), INDFUN( lcjac )
      REAL ( KIND = wp ) :: X( n ), C( m ), CJAC ( lcjac )

!  Compute the values of the constraint functions and their gradients
!  for constraints initially written in Standard Input Format (SIF).
!  The Jacobian must be stored in a sparse format.
! (Subroutine CCFG performs the same calculations for a dense Jacobian.)

!  CJAC  is an array which gives the values of the nonzeros of the
!        general constraint functions evaluated at X and V.
!        The i-th entry of CJAC gives the value of the derivative
!        with respect to variable INDVAR(i) of constraint function 
!        INDFUN(i) (i.e., INDFUN(i) = j > 0 indicates the j-th
!        general constraint function).

!  Based on the subroutines cfn.f and csgr.f by Nick Gould, which are
!  in turn based on the subroutine SBMIN by Conn, Gould and Toint.

!  Ingrid Bongartz
!  April 1992.

!  local variables.

      INTEGER :: i, j, iel, k, ig, ii, ig1, l, ll, icon, icnt
      INTEGER :: nin, nvarel, nelow, nelup, istrgv, iendgv
      INTEGER :: llo, llwrk, ifstat, igstat
      LOGICAL :: NONTRV
      EXTERNAL :: RANGE 
      REAL ( KIND = wp ) :: ftt, one, zero, gi, scalee
      PARAMETER ( zero = 0.0_wp, one = 1.0_wp )

      IF ( data%numcon == 0 ) RETURN

!  identify which elements are included in constraints. Use logical work 
!  vector to keep track of elements already included

      data%LOGIC( : data%nel ) = .FALSE.

!  now identify elements in first m constraint groups

      icnt = 0
      DO ig = 1, data%ng
        icon = data%KNDOFC( ig )
        IF ( icon > 0 .AND. icon <= m ) THEN
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

      DO ig = 1, data%ng
        ftt = zero

!  Consider only those groups in the constraints

        icon = data%KNDOFC( ig )
        IF ( icon > 0 .AND. icon <= m ) THEN
          ftt = - data%B( ig )

!  Include contributions from the linear element only if the variable belongs 
!  to the first n variables

          DO i = data%ISTADA( ig ), data%ISTADA( ig + 1 ) - 1
            j = data%ICNA( i )
            IF ( j <= n ) ftt = ftt + data%A( i ) * X( j )
          END DO

!  include the contributions from the nonlinear elements

          DO i = data%ISTADG( ig ), data%ISTADG( ig + 1 ) - 1
            ftt = ftt + data%ESCALE( i ) * data%FUVALS( data%IELING( i ) )
          END DO

!  record derivatives of trivial groups

          IF ( data%GXEQX( ig ) ) data%GVALS( ig, 2 ) = one 
        END IF
        data%FT( ig ) = ftt
      END DO

!  compute the group function values

!  all group functions are trivial

      IF ( data%altriv ) THEN
        data%GVALS( : data%ng, 1 ) = data%FT( : data%ng )
        data%GVALS( : data%ng, 2 ) = one
      ELSE

!  evaluate the group function values. Evaluate groups belonging to the first 
!  m constraints only

        icnt = 0
        DO ig = 1, data%ng
          icon = data%KNDOFC( ig )
          IF ( icon > 0 .AND. icon <= m ) THEN
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

!  compute the constraint function values

      DO ig = 1, data%ng
        i = data%KNDOFC( ig )
        IF ( i > 0 .AND. i <= m ) THEN
          IF ( data%GXEQX( ig ) ) THEN
            C( i ) = data%GSCALE( ig ) * data%FT( ig )
          ELSE 
            C( i ) = data%GSCALE( ig ) * data%GVALS( ig, 1 )
          END IF
        END IF
      END DO

!  evaluate the group derivative values

      IF ( grad ) THEN
        IF ( .NOT. data%altriv ) THEN
          CALL GROUP( data%GVALS, data%ng, data%FT, data%GPVALU, data%ng,      &
                      data%ITYPEG, data%ISTGP, data%ICALCF, data%ltypeg,       &
                      data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,       &
                      .TRUE., igstat )
          IF ( igstat /= 0 ) GO TO 930
        END IF

!  compute the gradient values.  Initialize the Jacobian as zero

         nnzj = 0 
         CJAC( : lcjac ) = zero

!  consider the ig-th group

         DO ig = 1, data%ng
           icon = data%KNDOFC( ig )

!  consider only those groups in the first m constraints

           IF ( icon == 0 .OR. icon > m ) CYCLE
           ig1 = ig + 1
           istrgv = data%IWORK( data%lstagv + ig )
           iendgv = data%IWORK( data%lstagv + ig1 ) - 1
           nelow = data%ISTADG( ig )
           nelup = data%ISTADG( ig1 ) - 1
           nontrv = .NOT. data%GXEQX( ig )

!  compute the first derivative of the group

           gi = data%GSCALE( ig )
           IF ( nontrv ) gi = gi  * data%GVALS( ig, 2 )
             data%WRK( data%IWORK( data%lsvgrp + istrgv :                      &
                                   data%lsvgrp + iendgv ) ) = zero

!  the group has nonlinear elements

           IF ( nelow <= nelup ) THEN

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
                 CALL RANGE( iel, .TRUE., data%FUVALS( k ),                    &
                             data%WRK( n + 1 ), nvarel, nin,                   &
                             data%ITYPEE( iel ), nin, nvarel )
!DIR$ IVDEP
                 DO i = 1, nvarel
                   j = data%IELVAR( l )
                   data%WRK( j ) = data%WRK( j ) + scalee * data%WRK( n + i )
                   l = l + 1
                 END DO
               ELSE

!  the iel-th element has no internal representation

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

               IF ( ll <= n ) THEN
                 nnzj = nnzj + 1
                 IF ( nnzj <= lcjac ) THEN
                   CJAC ( nnzj ) = gi * data%WRK( ll )
                   INDFUN( nnzj ) = icon
                   INDVAR( nnzj ) = ll
                 END IF
               END IF
             END DO

!  the group has only linear elements

           ELSE

!  include the contribution from the linear element

!DIR$ IVDEP
             DO k = data%ISTADA( ig ),data%ISTADA( ig1 ) - 1
               j = data%ICNA( k )
               data%WRK( j ) = data%WRK( j ) + data%A( k )
             END DO

!  allocate a gradient

!DIR$ IVDEP
             DO i = istrgv, iendgv
               ll = data%IWORK( data%lsvgrp + i )

!  include contributions from the first n variables only

               IF ( ll <= n ) THEN
                 nnzj = nnzj + 1
                 IF ( nnzj <= lcjac ) THEN
                   CJAC ( nnzj ) = gi * data%WRK( ll )
                   INDFUN( nnzj ) = icon
                   INDVAR( nnzj ) = ll
                 END IF
               END IF
             END DO
           END IF
         END DO

!  verify that the Jacobian can fit in the allotted space

         IF ( nnzj > lcjac ) THEN
           IF ( data%out > 0 ) WRITE( data%out, 2000 ) nnzj - lcjac 
           status = 2 ; RETURN
         END IF
      END IF

!  Update the counters for the report tool.

      data%nc2cf = data%nc2cf + data%pnc
      IF ( grad ) data%nc2cg = data%nc2cg + data%pnc
      status = 0
      RETURN

!  unsuccessful returns

  930 CONTINUE
      IF ( data%out > 0 ) WRITE( data%out,                                     &
        "( ' ** SUBROUTINE CCFSG: error flag raised during SIF evaluation' )" )
      status = 3
      RETURN

!  non-executable statements

 2000 FORMAT( /  ' ** SUBROUTINE CCFSG: array length lcjac too small.'         &
              /  ' -- Increase the parameter lcjac by at least ', I0,          &
                 ' and restart' )

!  end of subroutine CSCFG

      END SUBROUTINE CCFSG

