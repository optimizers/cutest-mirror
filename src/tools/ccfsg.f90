! THIS VERSION: CUTEST 1.0 - 28/11/2012 AT 08:20 GMT.

!-*-*-*-*-*-*-*-  C U T E S T    C C F S G    S U B R O U T I N E  -*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal authors: Ingrid Bongartz and Nick Gould

!  History -
!   fortran 77 version originally released in CUTE, April 1992
!   fortran 2003 version released in CUTEst, 28th November 2012

      SUBROUTINE CCFSG( data, status, n, m, X, C, nnzj, lj, J_val,             &
                        J_var, J_fun, grad )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      TYPE ( CUTEST_data_type ), INTENT( INOUT ) :: data
      INTEGER, INTENT( IN ) :: n, m, lj
      INTEGER, INTENT( OUT ) :: nnzj, status
      LOGICAL, INTENT( IN ) :: grad
      INTEGER, INTENT( OUT ), DIMENSION( lj ) :: J_var, J_fun
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: X
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( m ) :: C
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( lj ) :: J_val

!  ----------------------------------------------------------------------
!  Compute the values of the constraint functions and their gradients
!  for constraints initially written in Standard Input Format (SIF).
!  The Jacobian will be stored in a sparse co-ordinate format.
!  (Subroutine CCFG performs the same calculations for a dense Jacobian.)

!  J_val  is an array which gives the values of the nonzeros of the
!        general constraint functions evaluated at X and Y.
!        The i-th entry of J_val gives the value of the derivative
!        with respect to variable J_var(i) of constraint function 
!        J_fun(i) (i.e., J_fun(i) = j > 0 indicates the j-th
!        general constraint function).
!  ----------------------------------------------------------------------

!  local variables

      INTEGER :: i, j, iel, k, ig, ii, ig1, l, ll, icon, icnt
      INTEGER :: nin, nvarel, nelow, nelup, istrgv, iendgv, ifstat, igstat
      REAL ( KIND = wp ) :: ftt, gi, scalee
      EXTERNAL :: RANGE 

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
        ftt = 0.0_wp

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

          IF ( data%GXEQX( ig ) ) data%GVALS( ig, 2 ) = 1.0_wp
        END IF
        data%FT( ig ) = ftt
      END DO

!  compute the group function values

!  all group functions are trivial

      IF ( data%altriv ) THEN
        data%GVALS( : data%ng, 1 ) = data%FT( : data%ng )
        data%GVALS( : data%ng, 2 ) = 1.0_wp
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
         J_val( : lj ) = 0.0_wp

!  consider the ig-th group

         DO ig = 1, data%ng
           icon = data%KNDOFC( ig )

!  consider only those groups in the first m constraints

           IF ( icon == 0 .OR. icon > m ) CYCLE
           ig1 = ig + 1
           istrgv = data%ISTAGV( ig )
           iendgv = data%ISTAGV( ig1 ) - 1
           nelow = data%ISTADG( ig )
           nelup = data%ISTADG( ig1 ) - 1

!  compute the first derivative of the group

           gi = data%GSCALE( ig )
           IF ( .NOT. data%GXEQX( ig ) ) gi = gi  * data%GVALS( ig, 2 )
             data%W_ws( data%ISVGRP( istrgv : iendgv ) ) = 0.0_wp

!  the group has nonlinear elements

           IF ( nelow <= nelup ) THEN

!  loop over the group's nonlinear elements

             DO ii = nelow, nelup
               iel = data%IELING( ii )
               k = data%INTVAR( iel )
               l = data%ISTAEV( iel )
               nvarel = data%ISTAEV( iel + 1 ) - l
               scalee = data%ESCALE( ii )

!  the iel-th element has an internal representation

               IF ( data%INTREP( iel ) ) THEN
                 nin = data%INTVAR( iel + 1 ) - k
                 CALL RANGE( iel, .TRUE., data%FUVALS( k ), data%W_el,         &
                             nvarel, nin, data%ITYPEE( iel ), nin, nvarel )
!DIR$ IVDEP
                 DO i = 1, nvarel
                   j = data%IELVAR( l )
                   data%W_ws( j ) = data%W_ws( j ) + scalee * data%W_el( i )
                   l = l + 1
                 END DO

!  the iel-th element has no internal representation

               ELSE
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

!  include contributions from the first n variables only

               IF ( ll <= n ) THEN
                 nnzj = nnzj + 1
                 IF ( nnzj <= lj ) THEN
                   J_val ( nnzj ) = gi * data%W_ws( ll )
                   J_fun( nnzj ) = icon
                   J_var( nnzj ) = ll
                 END IF
               END IF
             END DO

!  the group has only linear elements

           ELSE

!  include the contribution from the linear element

!DIR$ IVDEP
             DO k = data%ISTADA( ig ),data%ISTADA( ig1 ) - 1
               j = data%ICNA( k )
               data%W_ws( j ) = data%W_ws( j ) + data%A( k )
             END DO

!  allocate a gradient

!DIR$ IVDEP
             DO i = istrgv, iendgv
               ll = data%ISVGRP( i )

!  include contributions from the first n variables only

               IF ( ll <= n ) THEN
                 nnzj = nnzj + 1
                 IF ( nnzj <= lj ) THEN
                   J_val ( nnzj ) = gi * data%W_ws( ll )
                   J_fun( nnzj ) = icon
                   J_var( nnzj ) = ll
                 END IF
               END IF
             END DO
           END IF
         END DO

!  verify that the Jacobian can fit in the allotted space

         IF ( nnzj > lj ) THEN
           IF ( data%out > 0 ) WRITE( data%out,                                &
             "( /,  ' ** SUBROUTINE CCFSG: array length lj too small', /,      &
            &   ' -- Increase the parameter lj to at least ', I0 )" ) nnzj
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

!  end of subroutine CSCFG

      END SUBROUTINE CCFSG

