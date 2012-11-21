! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CSGR ( data, status, n, m, grlagf, lv, V, X, nnzj, &
                        lcjac, CJAC, INDVAR, INDFUN )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n, m, lv, nnzj, lcjac
      INTEGER, INTENT( OUT ) :: status
      LOGICAL :: grlagf
      INTEGER :: INDVAR( LCJAC), INDFUN( lcjac )
      REAL ( KIND = wp ) :: X( n ), V ( lv ), CJAC ( lcjac )

!  Compute the gradients of the objective function and general
!  constraints of a function initially written in Standard
!  Input Format (SIF). The gradients are given in a sparse format.

!  CJAC	 is an array which gives the values of the nonzeros of the
!	 gradients of the objective, or Lagrangian, and general
!	 constraint functions evaluated  at X and V. The i-th entry
!	 of CJAC gives the value of the derivative with respect to
!	 variable INDVAR(i) of function INDFUN(i). INDFUN(i) = 0
!        indicates the objective function whenever grlagf is .FALSE.
!        or the Lagrangian function when grlagf is .TRUE., while
!        INDFUN(i) = j > 0 indicates the j-th general constraint
!        function.

!  Based on the minimization subroutine data%laNCELOT/SBMIN
!  by Conn, Gould and Toint.

!  Nick Gould, for CGT productions,
!  November 1991.

!  local variables.

      INTEGER :: i, j, iel, k, ig, ii, ig1, l, jj, ll, icon, ifstat, igstat
      INTEGER :: nin, nvarel, nelow, nelup, istrgv, iendgv, lgtemp
      REAL ( KIND = wp ) :: ftt, gi, scalee, gii
      LOGICAL :: nontrv
      EXTERNAL :: RANGE 

!  there are non-trivial group functions

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

!  evaluate the element function gradient values

      CALL ELFUN( data%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,          &
                  data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,          &
                  data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,           &
                  data%lelvar, data%lntvar, data%lstadh, data%lstep,           &
                  data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,          &
                  2, ifstat )

!  compute the group argument values ft

      DO ig = 1, data%ng
        ftt = - data%B( ig )

!  include the contribution from the linear element.

        DO j = data%ISTADA( ig ), data%ISTADA( ig + 1 ) - 1
          ftt = ftt + data%A( j ) * X( data%ICNA( j ) )
        END DO

!  include the contributions from the nonlinear elements.

        DO j = data%ISTADG( ig ), data%ISTADG( ig + 1 ) - 1
          ftt = ftt + data%ESCALE( j ) * data%FUVALS( data%IELING( j ) )
        END DO
        data%FT( ig ) = ftt

!  Record the derivatives of trivial groups.

        IF ( data%GXEQX( ig ) ) data%GVALS( ig, 2 ) = 1.0_wp
      END DO

!  evaluate the group derivative values.

      IF ( .NOT. data%altriv )                                                 &
        CALL GROUP( data%GVALS, data%ng, data%FT, data%GPVALU, data%ng,        &
                    data%ITYPEG, data%ISTGP, data%ICALCF, data%ltypeg,         &
                    data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,         &
                    .TRUE., igstat )

!  compute the gradient values. Initialize the gradient of the objective 
!  function as zero

      nnzj = 0
      IF ( data%numcon > 0 ) THEN
        lgtemp = WRK + n + data%maxsel
        DO j = 1, n
          data%WRK( lgtemp + j ) = 0.0_wp
        END DO

!  consider the IG-th group

        DO ig = 1, data%ng
          ig1 = ig + 1
          icon = data%KNDOFC( ig )
          istrgv = data%IWORK( data%lstagv + ig )
          iendgv = data%IWORK( data%lstagv + ig1 ) - 1
          nelow = data%ISTADG( ig )
          nelup = data%ISTADG( ig1 ) - 1
          nontrv = .NOT. data%GXEQX( ig )

!  compute the first derivative of the group

          gi = data%GSCALE( ig )
          IF ( icon == 0 ) THEN
            gii = gi
          ELSE
            IF ( grlagf ) gii = gi * V( data%KNDOFC( ig ) )
          END IF
          IF ( nontrv ) THEN
            gi = gi  * data%GVALS( ig, 2 )
            IF ( grlagf ) gii = gii * data%GVALS( ig, 2 )
          END IF
          data%WRK( data%IWORK( data%lsvgrp + istrgv :                         &
                                data%lsvgrp + iendgv ) ) = 0.0_wp

!  this is the first gradient evaluation or the group has nonlinear elements

          IF ( data%firstg .OR. nelow <= nelup ) THEN

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
                CALL RANGE ( iel, .TRUE., data%FUVALS( k ),                    &
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
                  k = k + 1
                  l = l + 1
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

!  the group belongs to the objective function

              IF ( icon == 0 ) THEN
                data%WRK( lgtemp + ll )                                        &
                  = data%WRK( lgtemp + ll ) + gi * data%WRK( ll )

!  the group defines a constraint

              ELSE
                nnzj = nnzj + 1
                IF ( nnzj <= lcjac ) THEN
                   CJAC ( nnzj ) = gi * data%WRK( ll )
                   INDFUN( nnzj ) = icon
                   INDVAR( nnzj ) = ll
                END IF
                IF ( grlagf ) data%WRK( lgtemp + ll )                          &
                  = data%WRK( lgtemp + ll ) + gii * data%WRK( ll )
              END IF

!  if the group is non-trivial, also store the nonzero entries of the
!  gradient of the function in GRJAC

              IF ( nontrv ) THEN
                jj = data%IWORK( data%lstajc + ll )
                data%FUVALS( data%lgrjac + jj ) = data%WRK( ll )

!  increment the address for the next nonzero in the column of the Jacobian 
!  for variable ll

                data%IWORK( data%lstajc + ll ) = jj + 1
              END IF
            END DO

!  this is not the first gradient evaluation and there is only a linear element

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

!  the group belongs to the objective function

              IF ( icon == 0 ) THEN
                data%WRK( lgtemp + ll )                                        &
                  = data%WRK( lgtemp + ll ) +  gi * data%WRK( ll )

!  the group defines a constraint

               ELSE
                 nnzj = nnzj + 1
                 IF ( nnzj <= lcjac ) THEN
                   CJAC ( nnzj ) = gi * data%WRK( ll )
                   INDFUN( nnzj ) = icon
                   INDVAR( nnzj ) = ll
                 END IF
                 IF ( grlagf ) data%WRK( lgtemp + ll )                         &
                   = data%WRK( lgtemp + ll ) + gii * data%WRK( ll )
               END IF

!  increment the address for the next nonzero in the column of the Jacobian 
!  for variable ll

               IF ( nontrv ) THEN
                 jj = data%IWORK( data%lstajc + ll )
                 data%IWORK( data%lstajc + ll ) = jj + 1
               END IF
            END DO
          END IF
        END DO

!  reset the starting addresses for the lists of groups using each variable to 
!  their values on entry

        DO i = n, 2, - 1
          data%IWORK( data%lstajc + i ) = data%IWORK( data%lstajc + i - 1 )
        END DO
        data%IWORK( data%lstajc + 1 ) = 1

!  transfer the gradient of the objective function to the sparse storage scheme

        DO i = 1, n
          nnzj = nnzj + 1
          IF ( nnzj <= lcjac ) THEN
            CJAC ( nnzj ) = data%WRK( lgtemp + i )
            INDFUN( nnzj ) = 0
            INDVAR( nnzj ) = i
          END IF
        END DO
      ELSE

!  compute the gradient value

        CALL ELGRD( n, data%ng, data%firstg, data%ICNA( 1 ), data%licna, &
                    data%ISTADA( 1 ), data%lstada, data%IELING( 1 ), &
                    data%leling, data%ISTADG( 1 ), data%lstadg, &
                    data%ITYPEE( 1 ), data%lintre, &
                    data%ISTAEV( 1 ), data%lstaev, data%IELVAR( 1 ), &
                    data%lelvar, data%INTVAR( 1 ), data%lntvar, &
                    data%IWORK( data%lsvgrp + 1 ), &
                    data%lnvgrp, data%IWORK( data%lstajc + 1 ), data%lnstjc, &
                    data%IWORK( data%lstagv + 1 ), data%lnstgv, &
                    data%A( 1 ), data%la, &
                    data%GVALS( : , 2 ), data%lgvals, &
                    data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ), &
                    data%GSCALE( 1 ), data%lgscal, &
                    data%ESCALE( 1 ), data%lescal, &
                    data%FUVALS( data%lgrjac + 1 ), &
                    data%lngrjc, data%WRK( 1 ), data%WRK( n + 1 ), &
                    data%maxsel, &
                    data%GXEQX( 1 ), data%lgxeqx, &
                    data%INTREP( 1 ), data%lintre, RANGE )

!  transfer the gradient of the objective function to the sparse storage scheme

        DO i = 1, n
          nnzj = nnzj + 1
          IF ( nnzj <= lcjac ) THEN
            CJAC ( nnzj ) = data%FUVALS( data%lggfx + i )
            INDFUN( nnzj ) = 0
            INDVAR( nnzj ) = i
          END IF
        END DO
      END IF
      data%firstg = .FALSE.

!  verify that the Jacobian can fit in the alloted space

      IF ( nnzj > lcjac ) THEN
        IF ( data%out > 0 ) WRITE( data%out, 2000 ) nnzj - lcjac 
        status = 2 ; RETURN
      END IF

!  update the counters for the report tool

      data%nc2og = data%nc2og + 1
      data%nc2cg = data%nc2cg + data%pnc
      status = 0
      RETURN

!  non-executable statements

 2000 FORMAT( /, ' ** SUBROUTINE CSGR: array length lcjac too small.',  &
              /, ' -- Increase the parameter lcjac by at least ', I0, &
                 ' and restart.' )

!  end of CSGR

      END
