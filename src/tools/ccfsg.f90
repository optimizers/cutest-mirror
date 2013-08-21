! THIS VERSION: CUTEST 1.1 - 22/08/2013 AT 12:55 GMT.

!-*-*-*-*-*-  C U T E S T  C I N T _ C C F S G    S U B R O U T I N E  -*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 2003 version released in CUTEst, 21st August 2013

      SUBROUTINE CUTEST_Cint_ccfsg( status, n, m, X, C,                        &
                               nnzj, lj, J_val, J_var, J_fun, grad )
      USE CUTEST
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_Bool
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      INTEGER, INTENT( IN ) :: n, m, lj
      INTEGER, INTENT( OUT ) :: nnzj, status
      LOGICAL ( KIND = C_Bool ), INTENT( IN ) :: grad
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

      LOGICAL :: grad_fortran

      grad_fortran = grad
      CALL CUTEST_ccfsg( status, n, m, X, C,                                   &
                         nnzj, lj, J_val, J_var, J_fun, grad_fortran )

      RETURN

!  end of subroutine CUTEST_Cint_cscfg

      END SUBROUTINE CUTEST_Cint_ccfsg

!-*-*-*-*-*-*-*-  C U T E S T    C C F S G    S U B R O U T I N E  -*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 2003 version released in CUTEst, 29th December 2012

      SUBROUTINE CUTEST_ccfsg( status, n, m, X, C,                             &
                               nnzj, lj, J_val, J_var, J_fun, grad )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

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

      CALL CUTEST_ccfsg_threadsafe( CUTEST_data_global,                        &
                                    CUTEST_work_global( 1 ),                   &
                                    status, n, m, X, C,                        &
                                    nnzj, lj, J_val, J_var, J_fun, grad )
      RETURN

!  end of subroutine CUTEST_cscfg

      END SUBROUTINE CUTEST_ccfsg

!-*-*-  C U T E S T    C C F S G _ t h r e a d e d   S U B R O U T I N E  -*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 2003 version released in CUTEst, 29th December 2012

      SUBROUTINE CUTEST_ccfsg_threaded( status, n, m, X, C, nnzj, lj,          &
                                        J_val, J_var, J_fun, grad, thread )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      INTEGER, INTENT( IN ) :: n, m, lj, thread
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

!  check that the specified thread is within range

      IF ( thread < 1 .OR. thread > CUTEST_data_global%threads ) THEN
        IF ( CUTEST_data_global%out > 0 )                                      &
          WRITE( CUTEST_data_global%out, "( ' ** CUTEST error: thread ', I0,   &
         &  ' out of range [1,', I0, ']' )" ) thread, CUTEST_data_global%threads
        status = 4 ; RETURN
      END IF

!  evaluate using specified thread

      CALL CUTEST_ccfsg_threadsafe( CUTEST_data_global,                        &
                                    CUTEST_work_global( thread ),              &
                                    status, n, m, X, C,                        &
                                    nnzj, lj, J_val, J_var, J_fun, grad )
      RETURN

!  end of subroutine CUTEST_cscfg_threaded

      END SUBROUTINE CUTEST_ccfsg_threaded

!-*-  C U T E S T    C C F S G _ t h r e a d s a f e   S U B R O U T I N E  -*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal authors: Ingrid Bongartz and Nick Gould

!  History -
!   fortran 77 version originally released in CUTE, April 1992
!   fortran 2003 version released in CUTEst, 28th November 2012

      SUBROUTINE CUTEST_ccfsg_threadsafe( data, work, status, n, m, X, C,      &
                                          nnzj, lj, J_val, J_var, J_fun, grad )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      TYPE ( CUTEST_data_type ), INTENT( IN ) :: data
      TYPE ( CUTEST_work_type ), INTENT( INOUT ) :: work
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

      work%LOGIC( : data%nel ) = .FALSE.

!  now identify elements in first m constraint groups

      icnt = 0
      DO ig = 1, data%ng
        icon = data%KNDOFC( ig )
        IF ( icon > 0 .AND. icon <= m ) THEN
          nelow = data%ISTADG( ig )
          nelup = data%ISTADG( ig + 1 ) - 1
          DO ii = nelow, nelup
            iel = data%IELING( ii )
            IF ( .NOT. work%LOGIC( iel ) ) THEN
              work%LOGIC( iel ) = .TRUE.
              icnt = icnt + 1
              work%ICALCF( icnt ) = iel
            END IF
          END DO
        END IF
      END DO

!  evaluate the element function values

      CALL ELFUN( work%FUVALS, X, data%EPVALU, icnt, data%ITYPEE,              &
                  data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,          &
                  data%ISTEP, work%ICALCF, data%ltypee, data%lstaev,           &
                  data%lelvar, data%lntvar, data%lstadh, data%lstep,           &
                  data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,          &
                  1, ifstat )
      IF ( ifstat /= 0 ) GO TO 930

!  evaluate the element function derivatives

      IF ( grad )                                                              &
        CALL ELFUN( work%FUVALS, X, data%EPVALU, icnt, data%ITYPEE,            &
                    data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,        &
                    data%ISTEP, work%ICALCF, data%ltypee, data%lstaev,         &
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
            ftt = ftt + data%ESCALE( i ) * work%FUVALS( data%IELING( i ) )
          END DO

!  record derivatives of trivial groups

          IF ( data%GXEQX( ig ) ) work%GVALS( ig, 2 ) = 1.0_wp
        END IF
        work%FT( ig ) = ftt
      END DO

!  compute the group function values

!  all group functions are trivial

      IF ( data%altriv ) THEN
        work%GVALS( : data%ng, 1 ) = work%FT( : data%ng )
        work%GVALS( : data%ng, 2 ) = 1.0_wp
      ELSE

!  evaluate the group function values. Evaluate groups belonging to the first 
!  m constraints only

        icnt = 0
        DO ig = 1, data%ng
          icon = data%KNDOFC( ig )
          IF ( icon > 0 .AND. icon <= m ) THEN
            icnt = icnt + 1
            work%ICALCF( icnt ) = ig
          END IF 
        END DO
        CALL GROUP( work%GVALS, data%ng, work%FT, data%GPVALU, icnt,           &
                    data%ITYPEG, data%ISTGP, work%ICALCF, data%ltypeg,         &
                    data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,         &
                    .FALSE., igstat )
        IF ( igstat /= 0 ) GO TO 930
      END IF

!  compute the constraint function values

      DO ig = 1, data%ng
        i = data%KNDOFC( ig )
        IF ( i > 0 .AND. i <= m ) THEN
          IF ( data%GXEQX( ig ) ) THEN
            C( i ) = data%GSCALE( ig ) * work%FT( ig )
          ELSE 
            C( i ) = data%GSCALE( ig ) * work%GVALS( ig, 1 )
          END IF
        END IF
      END DO

!  evaluate the group derivative values

      IF ( grad ) THEN
        IF ( .NOT. data%altriv ) THEN
          CALL GROUP( work%GVALS, data%ng, work%FT, data%GPVALU, icnt,         &
                      data%ITYPEG, data%ISTGP, work%ICALCF, data%ltypeg,       &
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
           IF ( .NOT. data%GXEQX( ig ) ) gi = gi  * work%GVALS( ig, 2 )
             work%W_ws( data%ISVGRP( istrgv : iendgv ) ) = 0.0_wp

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
                 CALL RANGE( iel, .TRUE., work%FUVALS( k ), work%W_el,         &
                             nvarel, nin, data%ITYPEE( iel ), nin, nvarel )
!DIR$ IVDEP
                 DO i = 1, nvarel
                   j = data%IELVAR( l )
                   work%W_ws( j ) = work%W_ws( j ) + scalee * work%W_el( i )
                   l = l + 1
                 END DO

!  the iel-th element has no internal representation

               ELSE
!DIR$ IVDEP
                 DO i = 1, nvarel
                   j = data%IELVAR( l )
                   work%W_ws( j ) = work%W_ws( j ) + scalee * work%FUVALS( k )
                   k = k + 1 ; l = l + 1
                 END DO
               END IF
             END DO

!  include the contribution from the linear element

!DIR$ IVDEP
             DO k = data%ISTADA( ig ), data%ISTADA( ig1 ) - 1
               j = data%ICNA( k )
               work%W_ws( j ) = work%W_ws( j ) + data%A( k )
             END DO

!  allocate a gradient

!DIR$ IVDEP
             DO i = istrgv, iendgv
               ll = data%ISVGRP( i )

!  include contributions from the first n variables only

               IF ( ll <= n ) THEN
                 nnzj = nnzj + 1
                 IF ( nnzj <= lj ) THEN
                   J_val ( nnzj ) = gi * work%W_ws( ll )
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
               work%W_ws( j ) = work%W_ws( j ) + data%A( k )
             END DO

!  allocate a gradient

!DIR$ IVDEP
             DO i = istrgv, iendgv
               ll = data%ISVGRP( i )

!  include contributions from the first n variables only

               IF ( ll <= n ) THEN
                 nnzj = nnzj + 1
                 IF ( nnzj <= lj ) THEN
                   J_val ( nnzj ) = gi * work%W_ws( ll )
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

      work%nc2cf = work%nc2cf + work%pnc
      IF ( grad ) work%nc2cg = work%nc2cg + work%pnc
      status = 0
      RETURN

!  unsuccessful returns

  930 CONTINUE
      IF ( data%out > 0 ) WRITE( data%out,                                     &
        "( ' ** SUBROUTINE CCFSG: error flag raised during SIF evaluation' )" )
      status = 3
      RETURN

!  end of subroutine CUTEST_cscfg_threadsafe

      END SUBROUTINE CUTEST_ccfsg_threadsafe

