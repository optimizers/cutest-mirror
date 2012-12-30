! THIS VERSION: CUTEST 1.0 - 29/12/2012 AT 11:30 GMT.

!-*-*-*-*-*-*-*-  C U T E S T    C S G R   S U B R O U T I N E  -*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 2003 version released in CUTEst, 29th December 2012

      SUBROUTINE CUTEST_csgr( status, n, m, X, Y, grlagf,                      &
                              nnzj, lj, J_val, J_var, J_fun )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      INTEGER, INTENT( IN ) :: n, m, lj
      INTEGER, INTENT( OUT ) :: nnzj, status
      LOGICAL, INTENT( IN ) :: grlagf
      INTEGER, INTENT( OUT ), DIMENSION( lj ) :: J_var, J_fun
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: X
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( m ) :: Y
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( lj ) :: J_val

!  ----------------------------------------------------------------
!  compute the gradients of the objective function and general
!  constraints of a function initially written in Standard
!  Input Format (SIF). The gradients are given in a sparse format

!  J_val is an array which gives the values of the nonzeros of the
!	 gradients of the objective, or Lagrangian, and general
!	 constraint functions evaluated  at X and Y. The i-th entry
!	 of J_val gives the value of the derivative with respect to
!	 variable J_var(i) of function J_fun(i). J_fun(i) = 0
!        indicates the objective function whenever grlagf is .FALSE.
!        or the Lagrangian function when grlagf is .TRUE., while
!        J_fun(i) = j > 0 indicates the j-th general constraint
!        function
!  ----------------------------------------------------------------

      CALL CUTEST_csgr_threadsafe( CUTEST_data_global,                         &
                                   CUTEST_work_global( 1 ),                    &
                                   status, n, m, X, Y, grlagf,                 &
                                   nnzj, lj, J_val, J_var, J_fun )
      RETURN

!  end of subroutine CUTEST_csgr

      END SUBROUTINE CUTEST_csgr

!-*-*-*-  C U T E S T    C S G R _ t h r e a d e d  S U B R O U T I N E  -*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 2003 version released in CUTEst, 29th December 2012

      SUBROUTINE CUTEST_csgr_threaded( status, n, m, X, Y, grlagf,             &
                                       nnzj, lj, J_val, J_var, J_fun, thread )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      INTEGER, INTENT( IN ) :: n, m, lj, thread
      INTEGER, INTENT( OUT ) :: nnzj, status
      LOGICAL, INTENT( IN ) :: grlagf
      INTEGER, INTENT( OUT ), DIMENSION( lj ) :: J_var, J_fun
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: X
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( m ) :: Y
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( lj ) :: J_val

!  ----------------------------------------------------------------
!  compute the gradients of the objective function and general
!  constraints of a function initially written in Standard
!  Input Format (SIF). The gradients are given in a sparse format

!  J_val is an array which gives the values of the nonzeros of the
!	 gradients of the objective, or Lagrangian, and general
!	 constraint functions evaluated  at X and Y. The i-th entry
!	 of J_val gives the value of the derivative with respect to
!	 variable J_var(i) of function J_fun(i). J_fun(i) = 0
!        indicates the objective function whenever grlagf is .FALSE.
!        or the Lagrangian function when grlagf is .TRUE., while
!        J_fun(i) = j > 0 indicates the j-th general constraint
!        function
!  ----------------------------------------------------------------

      CALL CUTEST_csgr_threadsafe( CUTEST_data_global,                         &
                                   CUTEST_work_global( thread ),               &
                                   status, n, m, X, Y, grlagf,                 &
                                   nnzj, lj, J_val, J_var, J_fun )
      RETURN

!  end of subroutine CUTEST_csgr_threaded

      END SUBROUTINE CUTEST_csgr_threaded

!-*-*-  C U T E S T    C S G R _ t h r e a d s a f e  S U B R O U T I N E  -*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 77 version originally released in CUTE, November 1991
!   fortran 2003 version released in CUTEst, 28th November 2012

      SUBROUTINE CUTEST_csgr_threadsafe( data, work, status, n, m, X, Y,       &
                                         grlagf, nnzj, lj, J_val, J_var, J_fun )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      TYPE ( CUTEST_data_type ), INTENT( IN ) :: data
      TYPE ( CUTEST_work_type ), INTENT( INOUT ) :: work
      INTEGER, INTENT( IN ) :: n, m, lj
      INTEGER, INTENT( OUT ) :: nnzj, status
      LOGICAL, INTENT( IN ) :: grlagf
      INTEGER, INTENT( OUT ), DIMENSION( lj ) :: J_var, J_fun
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: X
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( m ) :: Y
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( lj ) :: J_val

!  ----------------------------------------------------------------
!  compute the gradients of the objective function and general
!  constraints of a function initially written in Standard
!  Input Format (SIF). The gradients are given in a sparse format

!  J_val is an array which gives the values of the nonzeros of the
!	 gradients of the objective, or Lagrangian, and general
!	 constraint functions evaluated  at X and Y. The i-th entry
!	 of J_val gives the value of the derivative with respect to
!	 variable J_var(i) of function J_fun(i). J_fun(i) = 0
!        indicates the objective function whenever grlagf is .FALSE.
!        or the Lagrangian function when grlagf is .TRUE., while
!        J_fun(i) = j > 0 indicates the j-th general constraint
!        function
!  ----------------------------------------------------------------

!  local variables

      INTEGER :: i, j, iel, k, ig, ii, ig1, l, jj, ll, icon, ifstat, igstat
      INTEGER :: nin, nvarel, nelow, nelup, istrgv, iendgv
      REAL ( KIND = wp ) :: ftt, gi, scalee, gii
      LOGICAL :: nontrv
      EXTERNAL :: RANGE 

!  there are non-trivial group functions

      DO i = 1, MAX( data%nel, data%ng )
        work%ICALCF( i ) = i
      END DO

!  evaluate the element function values

      CALL ELFUN( work%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,          &
                  data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,          &
                  data%ISTEP, work%ICALCF, data%ltypee, data%lstaev,           &
                  data%lelvar, data%lntvar, data%lstadh, data%lstep,           &
                  data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,          &
                  1, ifstat )
      IF ( ifstat /= 0 ) GO TO 930

!  evaluate the element function gradient values

      CALL ELFUN( work%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,          &
                  data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,          &
                  data%ISTEP, work%ICALCF, data%ltypee, data%lstaev,           &
                  data%lelvar, data%lntvar, data%lstadh, data%lstep,           &
                  data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,          &
                  2, ifstat )
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
          ftt = ftt + data%ESCALE( j ) * work%FUVALS( data%IELING( j ) )
        END DO
        work%FT( ig ) = ftt

!  record the derivatives of trivial groups

        IF ( data%GXEQX( ig ) ) work%GVALS( ig, 2 ) = 1.0_wp
      END DO

!  evaluate the group derivative values

      IF ( .NOT. data%altriv ) THEN
        CALL GROUP( work%GVALS, data%ng, work%FT, data%GPVALU, data%ng,        &
                    data%ITYPEG, data%ISTGP, work%ICALCF, data%ltypeg,         &
                    data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,         &
                    .TRUE., igstat )
        IF ( igstat /= 0 ) GO TO 930
      END IF

!  compute the gradient values. Initialize the gradient of the objective 
!  function as zero

      nnzj = 0
      IF ( data%numcon > 0 ) THEN
        work%G_temp( : n ) = 0.0_wp

!  consider the ig-th group

        DO ig = 1, data%ng
          ig1 = ig + 1
          icon = data%KNDOFC( ig )
          istrgv = data%ISTAGV( ig )
          iendgv = data%ISTAGV( ig1 ) - 1
          nelow = data%ISTADG( ig )
          nelup = data%ISTADG( ig1 ) - 1
          nontrv = .NOT. data%GXEQX( ig )

!  compute the first derivative of the group

          gi = data%GSCALE( ig )
          IF ( icon == 0 ) THEN
            gii = gi
          ELSE
            IF ( grlagf ) gii = gi * Y( data%KNDOFC( ig ) )
          END IF
          IF ( nontrv ) THEN
            gi = gi * work%GVALS( ig, 2 )
            IF ( grlagf ) gii = gii * work%GVALS( ig, 2 )
          END IF
          work%W_ws( data%ISVGRP( istrgv : iendgv ) ) = 0.0_wp

!  this is the first gradient evaluation or the group has nonlinear elements

          IF ( work%firstg .OR. nelow <= nelup ) THEN

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
                CALL RANGE( iel, .TRUE., work%FUVALS( k ), work%W_el,          &
                            nvarel, nin, data%ITYPEE( iel ), nin, nvarel )
!DIR$ IVDEP
                DO i = 1, nvarel
                  j = data%IELVAR( l )
                  work%W_ws( j ) = work%W_ws( j ) + scalee * work%W_el( i )
                  l = l + 1
                END DO
              ELSE

!  the iel-th element has no internal representation

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

!  the group belongs to the objective function

              IF ( icon == 0 ) THEN
                work%G_temp( ll ) = work%G_temp( ll ) + gi * work%W_ws( ll )

!  the group defines a constraint

              ELSE
                nnzj = nnzj + 1
                IF ( nnzj <= lj ) THEN
                   J_val ( nnzj ) = gi * work%W_ws( ll )
                   J_fun( nnzj ) = icon
                   J_var( nnzj ) = ll
                END IF
                IF ( grlagf )                                                  &
                  work%G_temp( ll ) = work%G_temp( ll ) + gii * work%W_ws( ll )
              END IF

!  if the group is non-trivial, also store the nonzero entries of the
!  gradient of the function in GRJAC

              IF ( nontrv ) THEN
                jj = work%ISTAJC( ll )
                work%FUVALS( data%lgrjac + jj ) = work%W_ws( ll )

!  increment the address for the next nonzero in the column of the Jacobian 
!  for variable ll

                work%ISTAJC( ll ) = jj + 1
              END IF
            END DO

!  this is not the first gradient evaluation and there is only a linear element

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

!  the group belongs to the objective function

              IF ( icon == 0 ) THEN
                work%G_temp( ll ) = work%G_temp( ll ) +  gi * work%W_ws( ll )

!  the group defines a constraint

               ELSE
                 nnzj = nnzj + 1
                 IF ( nnzj <= lj ) THEN
                   J_val ( nnzj ) = gi * work%W_ws( ll )
                   J_fun( nnzj ) = icon
                   J_var( nnzj ) = ll
                 END IF
                 IF ( grlagf )                                                 &
                   work%G_temp( ll ) = work%G_temp( ll ) + gii * work%W_ws( ll )
               END IF

!  increment the address for the next nonzero in the column of the Jacobian 
!  for variable ll

               IF ( nontrv ) THEN
                 jj = work%ISTAJC( ll )
                 work%ISTAJC( ll ) = jj + 1
               END IF
            END DO
          END IF
        END DO

!  reset the starting addresses for the lists of groups using each variable to 
!  their values on entry

        DO i = n, 2, - 1
          work%ISTAJC( i ) = work%ISTAJC( i - 1 )
        END DO
        work%ISTAJC( 1 ) = 1

!  transfer the gradient of the objective function to the sparse storage scheme

        DO i = 1, n
          nnzj = nnzj + 1
          IF ( nnzj <= lj ) THEN
            J_val ( nnzj ) = work%G_temp( i )
            J_fun( nnzj ) = 0
            J_var( nnzj ) = i
          END IF
        END DO

!  compute the gradient value

      ELSE
        CALL CUTEST_form_gradients( n, data%ng, data%nel, data%ntotel,         &
               data%nvrels, data%nnza, data%nvargp, work%firstg, data%ICNA,    &
               data%ISTADA, data%IELING, data%ISTADG, data%ISTAEV,             &
               data%IELVAR, data%INTVAR, data%A, work%GVALS( : , 2 ),          &
               work%FUVALS, data%lnguvl, work%FUVALS( data%lggfx + 1 ),        &
               data%GSCALE, data%ESCALE, work%FUVALS( data%lgrjac + 1 ),       &
               data%GXEQX, data%INTREP, data%ISVGRP, data%ISTAGV, data%ITYPEE, &
               work%ISTAJC, work%W_ws, work%W_el, RANGE )

!  transfer the gradient of the objective function to the sparse storage scheme

        DO i = 1, n
          nnzj = nnzj + 1
          IF ( nnzj <= lj ) THEN
            J_val ( nnzj ) = work%FUVALS( data%lggfx + i )
            J_fun( nnzj ) = 0
            J_var( nnzj ) = i
          END IF
        END DO
      END IF
      work%firstg = .FALSE.

!  verify that the Jacobian can fit in the alloted space

      IF ( nnzj > lj ) THEN
        IF ( data%out > 0 ) WRITE( data%out,                                   &
          "( /, ' ** SUBROUTINE CSGR: array length lj too small.',             &
         &    /, ' -- Increase the parameter lj to at least ', I0 )" ) nnzj
        status = 2 ; RETURN
      END IF

!  update the counters for the report tool

      work%nc2og = work%nc2og + 1
      work%nc2cg = work%nc2cg + work%pnc
      status = 0
      RETURN

!  unsuccessful returns

  930 CONTINUE
      IF ( data%out > 0 ) WRITE( data%out,                                     &
        "( ' ** SUBROUTINE CSGR: error flag raised during SIF evaluation' )" )
      status = 3
      RETURN

!  end of subroutine CUTEST_csgr_threadsafe

      END SUBROUTINE CUTEST_csgr_threadsafe
