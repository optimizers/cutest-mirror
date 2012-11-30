! THIS VERSION: CUTEST 1.0 - 27/11/2012 AT 16:15 GMT.

!-*-*-*-*-*-*-  C U T E S T    C S G R E H    S U B R O U T I N E  -*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 77 version originally released in CUTEr, November 1994
!   fortran 2003 version released in CUTEst, 27th November 2012

      SUBROUTINE CSGREH( data, status, n, m, X, Y, grlagf, nnzj, lj,           &
                         J_val, J_var, J_fun, ne, lhe_ptr, HE_row_ptr,         &
                         HE_val_ptr, lhe_row, HE_row, lhe_val, HE_val, byrows )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      TYPE ( CUTEST_data_type ), INTENT( INOUT ) :: data
      INTEGER, INTENT( IN ) :: n, m, lj, lhe_ptr, lhe_row, lhe_val
      INTEGER, INTENT( OUT ) :: ne, nnzj, status
      LOGICAL, INTENT( IN ) :: grlagf, byrows
      INTEGER, INTENT( OUT ), DIMENSION( lj ) :: J_var, J_fun
      INTEGER, INTENT( OUT ), DIMENSION( lhe_ptr ) :: HE_row_ptr, HE_val_ptr
      INTEGER, INTENT( OUT ), DIMENSION( lhe_row ) :: HE_row
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: X
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( m ) :: Y
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( lj ) :: J_val
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( lhe_val ) :: HE_val

!  ----------------------------------------------------------------------------
!  compute the constraint Jacobian in co-ordinate format and Hessian matrix 
!  of the Lagrangian function of a problem initially written in Standard 
!  Input Format (SIF)

!  the Hessian matrix is represented in "finite element format", i.e., 

!           ne
!      H = sum H_e, 
!          e=1

!  where each element H_i involves a small subset of the rows of H. H is stored
!  as a list of the row indices involved in each element and the upper triangle
!  of H_e (stored by rows or columns). Specifically,

!  ne (integer) number of elements
!  HE_row (integer array) a list of the row indices involved which each
!          element. Those for element e directly proceed those for 
!          element e + 1, e = 1, ..., ne-1
!  HE_row_ptr (integer array) pointers to the position in HE_row of the first 
!          row index in each element. HE_row_ptr(ne+1) points to the first 
!          empty location in IRPNHI
!  HE_val (real array) a list of the nonzeros in the upper triangle of
!          H_e, stored by rows, or by columns, for each element. Those 
!          for element i directly proceed those for element, e + 1, 
!          e = 1, ..., ne-1
!  HE_val_ptr (integer array) pointers to the position in HE_val of the first 
!          nonzero in each element. HE_val_ptr(ne+1) points to the first 
!          empty location in HE_val
!  byrows (logical) must be set .TRUE. if the upper triangle of each H_e is
!          to be stored by rows, and .FALSE. if it is to be stored by columns
!  ----------------------------------------------------------------------------

!  Local variables

      INTEGER :: i, j, iel, k, ig, ii, ig1, l, jj, ll, icon
      INTEGER :: nin, nvarel, nelow, nelup, istrgv, iendgv
      INTEGER :: lhe_row_int, lhe_val_int, ifstat, igstat, alloc_status
      LOGICAL :: nontrv
      REAL ( KIND = wp ) :: ftt, gi, scalee, gii
      CHARACTER ( LEN = 80 ) :: bad_alloc = REPEAT( ' ', 80 )
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
      IF ( ifstat /= 0 ) GO TO 930

!  evaluate the element function gradients and Hessians

      CALL ELFUN( data%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,          &
                  data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,          &
                  data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,           &
                  data%lelvar, data%lntvar, data%lstadh, data%lstep,           &
                  data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,          &
                  3, ifstat )
      IF ( ifstat /= 0 ) GO TO 930

!  compute the group argument values ft

      DO ig = 1, data%ng
        ftt = - data%B( ig )

!  include the contribution from the linear element

        DO j = data%ISTADA( ig ), data%ISTADA( ig + 1 ) - 1
          ftt = ftt + data%A( j ) * X( data%ICNA( j ) )
        END DO

!  include the contributions from the nonlinear elements.

        DO j = data%ISTADG( ig ), data%ISTADG( ig + 1 ) - 1
          ftt = ftt + data%ESCALE( j ) * data%FUVALS( data%IELING( j ) )
        END DO
        data%FT( ig ) = ftt

!  Record the derivatives of trivial groups.

        IF ( data%GXEQX( ig ) ) THEN
          data%GVALS( ig, 2 ) = 1.0_wp
          data%GVALS( ig, 3 ) = 0.0_wp
        END IF
      END DO

!  evaluate the group derivative values

      IF ( .NOT. data%altriv ) THEN
        CALL GROUP( data%GVALS, data%ng, data%FT, data%GPVALU, data%ng,        &
                    data%ITYPEG, data%ISTGP, data%ICALCF, data%ltypeg,         &
                    data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,         &
                    .TRUE., igstat )
        IF ( igstat /= 0 ) GO TO 930
      END IF

!  define the real work space needed for ELGRD. Ensure that there is 
!  sufficient space

!     IF ( data%lwk2 < data%ng ) THEN
!       IF ( data%out > 0 ) WRITE( data%out, 2000 )
!       status = 2 ; RETURN
!     END IF

!  change the group weightings to include the contributions from the
!  Lagrange multipliers

      IF ( data%numcon > 0 ) THEN
         DO ig = 1, data%ng
           i = data%KNDOFC( ig )
           IF ( i == 0 ) THEN
             data%GSCALE_used( ig ) = data%GSCALE( ig )
           ELSE
             data%GSCALE_used( ig ) = data%GSCALE( ig ) * Y( i )
           END IF
         END DO

!  compute the gradient values. Initialize the gradient of the objective 
!  function as zero

         nnzj = 0
         data%G_temp( : n ) = 0.0_wp

!  consider the IG-th group

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
           gii = data%GSCALE_used( ig )
           IF ( nontrv ) THEN
             gi = gi  * data%GVALS( ig, 2 )
             gii = gii * data%GVALS( ig, 2 )
           END IF
           data%W_ws( data%ISVGRP( istrgv : iendgv ) ) = 0.0_wp

!  This is the first gradient evaluation or the group has nonlinear
!  elements.

           IF ( data%firstg .OR. nelow <= nelup ) THEN

!  loop over the group's nonlinear elements

             DO ii = nelow, nelup
               iel = data%IELING( ii )
               k = data%INTVAR( iel )
               l = data%ISTAEV( iel )
               nvarel = data%ISTAEV( iel + 1 ) - l
               scalee = data%ESCALE( ii )
               IF ( data%INTREP( iel ) ) THEN

!  the IEL-th element has an internal representation

                 nin = data%INTVAR( iel + 1 ) - k
                 CALL RANGE( iel, .TRUE., data%FUVALS( k ), data%W_el,         &
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

!  the group belongs to the objective function

               IF ( icon == 0 ) THEN
                 data%G_temp( ll ) = data%G_temp( ll ) + gi * data%W_ws( ll )

!  the group defines a constraint

               ELSE
                 nnzj = nnzj + 1
                 IF ( nnzj <= lj ) THEN
                   J_val ( nnzj ) = gi * data%W_ws( ll )
                   J_fun( nnzj ) = icon
                   J_var( nnzj ) = ll
                 END IF
                 IF ( grlagf )                                                 &
                   data%G_temp( ll ) = data%G_temp( ll ) + gii * data%W_ws( ll )
               END IF

!  if the group is non-trivial, also store the nonzero entries of the
!  gradient of the function in GRJAC.

               IF ( nontrv ) THEN
                 jj = data%ISTAJC( ll )
                 data%FUVALS( data%lgrjac + jj ) = data%W_ws( ll )

!  increment the address for the next nonzero in the column of
!  the jacobian for variable ll

                 data%ISTAJC( ll ) = jj + 1
               END IF
             END DO

!  this is not the first gradient evaluation and there is only a linear element

           ELSE

!  include the contribution from the linear element

!DIR$ IVDEP
             DO k = data%ISTADA( ig ),data%ISTADA( ig1 ) - 1
               j = data%ICNA( k )
               data%W_ws( j ) = data%W_ws( j ) + data%A( k )
             END DO

!  Allocate a gradient.

!DIR$ IVDEP
             DO i = istrgv, iendgv
               ll = data%ISVGRP( i )

!  the group belongs to the objective function

               IF ( icon == 0 ) THEN
                 data%G_temp( ll ) = data%G_temp( ll ) + gi * data%W_ws( ll )

!  the group defines a constraint

               ELSE
                 nnzj = nnzj + 1
                 IF ( nnzj <= lj ) THEN
                   J_val ( nnzj ) = gi * data%W_ws( ll )
                   J_fun( nnzj ) = icon
                   J_var( nnzj ) = ll
                 END IF
                 IF ( grlagf )                                                &
                   data%G_temp( ll ) = data%G_temp( ll ) + gii * data%W_ws( ll )
               END IF

!  increment the address for the next nonzero in the column of the Jacobian 
!  for variable ll

               IF ( nontrv ) THEN
                 jj = data%ISTAJC( ll )
                 data%ISTAJC( ll ) = jj + 1
               END IF
             END DO
           END IF
         END DO

!  reset the starting addresses for the lists of groups using each variable to
!  their values on entry

         DO i = n, 2, - 1
           data%ISTAJC( i ) = data%ISTAJC( i - 1 )
         END DO
         data%ISTAJC( 1 ) = 1

!  transfer the gradient of the objective function to the sparse storage scheme

        DO i = 1, n
          nnzj = nnzj + 1
          IF ( nnzj <= lj ) THEN
            J_val ( nnzj ) = data%G_temp( i )
            J_fun( nnzj ) = 0
            J_var( nnzj ) = i
          END IF
        END DO
      ELSE

!  compute the gradient value

        CALL CUTEST_form_gradients( n, data%ng, data%nel, data%ntotel,         &
               data%nvrels, data%nnza, data%nvargp, data%firstg, data%ICNA,    &
               data%ISTADA, data%IELING, data%ISTADG, data%ISTAEV,             &
               data%IELVAR, data%INTVAR, data%A, data%GVALS( : , 2 ),          &
               data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ),        &
               data%GSCALE, data%ESCALE, data%FUVALS( data%lgrjac + 1 ),       &
               data%GXEQX, data%INTREP, data%ISVGRP, data%ISTAGV, data%ITYPEE, &
               data%ISTAJC, data%W_ws, data%W_el, RANGE )

!        CALL ELGRD( n, data%ng, data%firstg, data%ICNA( 1 ), data%licna, &
!                      data%ISTADA( 1 ), data%lstada, data%IELING( 1 ), &
!                      data%leling, data%ISTADG( 1 ), data%lstadg, &
!                      data%ITYPEE( 1 ), data%lintre, &
!                      data%ISTAEV( 1 ), data%lstaev, data%IELVAR( 1 ), &
!                      data%lelvar, data%INTVAR( 1 ), data%lntvar, &
!                      data%ISVGRP( 1 ), &
!                      data%lnvgrp, data%ISTAJC( 1 ), data%lnstjc,&
!                      data%ISTAGV( 1 ), data%lnstgv, &
!                      data%A( 1 ), data%la, &
!                      data%GVALS( : , 2 ), data%lgvals, &
!                      data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ),&
!                      data%GSCALE( 1 ), data%lgscal, &
!                      data%ESCALE( 1 ), data%lescal, &
!                      data%FUVALS( data%lgrjac + 1 ), &
!                      data%lngrjc, data%W_ws( 1 ), data%W_el( 1 ), &
!                      data%maxsel, &
!                      data%GXEQX( 1 ), data%lgxeqx, &
!                      data%INTREP( 1 ), data%lintre, RANGE )

!  transfer the gradient of the objective function to the sparse storage scheme

         nnzj = 0
         DO i = 1, n
           nnzj = nnzj + 1
           IF ( nnzj <= lj ) THEN
             J_val ( nnzj ) = data%FUVALS( data%lggfx + i )
             J_fun( nnzj ) = 0
             J_var( nnzj ) = i
           END IF
         END DO
      END IF
      data%firstg = .FALSE.

!  verify that the Jacobian can fit in the alloted space

      IF ( nnzj > lj ) THEN
        IF ( data%out > 0 ) WRITE( data%out,                                   &
          "( /, ' ** SUBROUTINE CSGREH: array length lj too small.',           &
         &    /, ' -- Increase the parameter lj to at least ', I0 )" ) nnzj
        status = 2 ; RETURN
      END IF

!  define the real work space needed for ASMBE. Ensure that there is 
!  sufficient space

!     IF ( data%numcon > 0 ) THEN
!        IF ( data%lwk2 < n + 3 * data%maxsel + data%ng ) THEN
!           IF ( data%out > 0 ) WRITE( data%out, 2000 )
!           status = 2 ; RETURN
!        END IF
!     ELSE
!        IF ( data%lwk2 < n + 3 * data%maxsel ) THEN
!           IF ( data%out > 0 ) WRITE( data%out, 2000 )
!           status = 2 ; RETURN
!        END IF
 !    END IF

!  define the integer work space needed for ASMBE.  Ensure that there is 
!  sufficient space

!     liwkh = data%liwk2 - n

!  assemble the Hessian

      lhe_row_int = lhe_row ; lhe_val_int = lhe_val
      IF ( data%numcon > 0 ) THEN
        CALL CUTEST_assemble_element_hessian(                                  &
                        n, data%ng, data%nel,data% ntotel, data%nvrels,        &
                        data%nnza, data%maxsel, data%nvargp,                   &
                        data%lnguvl, data%lnhuvl, data%ISTADH, data%ICNA,      &
                        data%ISTADA, data%INTVAR, data%IELVAR,                 &
                        data%IELING, data%ISTADG, data%ISTAEV,                 &
                        data%ISTAGV, data%ISVGRP, data%ITYPEE,                 &
                        data%A, data%FUVALS, data%FUVALS,                      &
                        data%GVALS( : , 2 ), data%GVALS( : , 3 ),              &
                        data%GSCALE_used, data%ESCALE,                         &
                        data%GXEQX, data%INTREP,                               &
                        data%IW_asmbl, data%W_ws, data%W_el, data%W_in,        &
                        data%H_el, data%H_in, RANGE, ne, lhe_row_int,          &
                        lhe_val_int, data%H_row, HE_row_ptr, data%H_val,       &
                        HE_val_ptr, byrows, 0, data%out, data%out,             &
                        data%io_buffer, alloc_status, bad_alloc, status )
      ELSE
        CALL CUTEST_assemble_element_hessian(                                  &
                        n, data%ng, data%nel,data% ntotel, data%nvrels,        &
                        data%nnza, data%maxsel, data%nvargp,                   &
                        data%lnguvl, data%lnhuvl, data%ISTADH, data%ICNA,      &
                        data%ISTADA, data%INTVAR, data%IELVAR,                 &
                        data%IELING, data%ISTADG, data%ISTAEV,                 &
                        data%ISTAGV, data%ISVGRP, data%ITYPEE,                 &
                        data%A, data%FUVALS, data%FUVALS,                      &
                        data%GVALS( : , 2 ), data%GVALS( : , 3 ),              &
                        data%GSCALE, data%ESCALE, data%GXEQX, data%INTREP,     &
                        data%IW_asmbl, data%W_ws, data%W_el, data%W_in,        &
                        data%H_el, data%H_in, RANGE, ne, lhe_row_int,          &
                        lhe_val_int, data%H_row, HE_row_ptr, data%H_val,       &
                        HE_val_ptr, byrows, 0, data%out, data%out,             &
                        data%io_buffer, alloc_status, bad_alloc, status )
      END IF 

!  check for errors in the assembly

      IF ( status > 0 ) RETURN

!  check that HE_row and HE_val are large enough

      IF ( lhe_row < HE_row_ptr( ne + 1 ) - 1 ) THEN
        IF ( data%out > 0 ) WRITE( data%out, "( ' ** SUBROUTINE UEH: ',        &
       &  'Increase the dimension of HE_row to ',  I0 )" )                     &
             HE_row_ptr( ne + 1 ) - 1 
        status = 2 ; RETURN
      END IF

      IF ( lhe_val < HE_val_ptr( ne + 1 ) - 1 ) THEN
        IF ( data%out > 0 ) WRITE( data%out, "( ' ** SUBROUTINE UEH: ',        &
       &  'Increase the dimension of HE_val to ',  I0 )" )                     &
             HE_val_ptr( ne + 1 ) - 1 
        status = 2 ; RETURN
      END IF

!  record the element Hessian

      HE_row( : HE_row_ptr( ne + 1 ) - 1 )                                     &
         = data%H_row( : HE_row_ptr( ne + 1 ) - 1 )
      HE_val( : HE_val_ptr( ne + 1 ) - 1 )                                     &
         = data%H_val( : HE_val_ptr( ne + 1 ) - 1 )

!  assemble the Hessian

!      IF ( data%numcon > 0 ) THEN
!      CALL ASMBE( n, data%ng, data%maxsel,  &
!                      data%ISTADH( 1 ), data%lstadh, &
!                      data%ICNA( 1 ), data%licna, &
!                      data%ISTADA( 1 ), data%lstada, &
!                      data%INTVAR( 1 ), data%lntvar, &
!                      data%IELVAR( 1 ), data%lelvar, &
!                      data%IELING( 1 ), data%leling, &
!                      data%ISTADG( 1 ), data%lstadg, &
!                      data%ISTAEV( 1 ), data%lstaev, &
!                      data%ISTAGV( 1 ), data%lnstgv, &
!                      data%ISVGRP( 1 ), data%lnvgrp, &
!                      data%IWORK( liwkh + 1 ), data%liwk2 - liwkh, &
!                      data%A( 1 ), data%la, data%FUVALS, &
!                      data%lnguvl, data%FUVALS, data%lnhuvl, &
!                      data%GVALS( : , 2 ), data%GVALS( : , 3 ), &
!                      data%W_ws( 1 ), data%ESCALE( 1 ), data%lescal, &
!                      data%W_ws( data%ng + 1 ), data%lwk2 - data%ng, &
!                      data%GXEQX( 1 ), data%lgxeqx, data%INTREP( 1 ), &
!                      data%lintre, data%ITYPEE( 1 ), data%lintre, RANGE, ne,  &
!                      HE_row, lhe_row, HE_row_ptr, &
!                      HE_val, lhe_val, HE_val_ptr, &
!                      BYROWS, 1, data%out, inform )
!      ELSE
!      CALL ASMBE( n, data%ng, data%maxsel,  &
!                      data%ISTADH( 1 ), data%lstadh, &
!                      data%ICNA( 1 ), data%licna, &
!                      data%ISTADA( 1 ), data%lstada, &
!                      data%INTVAR( 1 ), data%lntvar, &
!                      data%IELVAR( 1 ), data%lelvar, &
!                      data%IELING( 1 ), data%leling, &
!                      data%ISTADG( 1 ), data%lstadg, &
!                      data%ISTAEV( 1 ), data%lstaev, &
!                      data%ISTAGV( 1 ), data%lnstgv, &
!                      data%ISVGRP( 1 ), data%lnvgrp, &
!                      data%IWORK( liwkh + 1 ), data%liwk2 - liwkh, &
!                      data%A( 1 ), data%la, data%FUVALS, &
!                      data%lnguvl, data%FUVALS, data%lnhuvl, &
!                      data%GVALS( : , 2 ), data%GVALS( : , 3 ), &
!                      data%GSCALE( 1 ), data%ESCALE( 1 ), data%lescal, &
!                      data%W_ws( 1 ), data%lwk2 - data%ng, &
!                      data%GXEQX( 1 ), data%lgxeqx, data%INTREP( 1 ), &
!                      data%lintre, data%ITYPEE( 1 ), data%lintre, RANGE, ne,  &
!                      HE_row, lhe_row, HE_row_ptr, &
!                      HE_val, lhe_val, HE_val_ptr, &
!                      BYROWS, 1, data%out, inform )
!      END IF

!  check that there is room for the elements

!     IF ( inform > 0 ) THEN
!       IF ( data%out > 0 ) WRITE( data%out, 2020 )
!       status = 2 ; RETURN
!     END IF

!  update the counters for the report tool

      data%nc2cg = data%nc2og + 1
      data%nc2oh = data%nc2oh + 1
      data%nc2cg = data%nc2cg + data%pnc
      data%nc2ch = data%nc2ch + data%pnc
      status = 0
      RETURN

!  unsuccessful returns

  930 CONTINUE
      IF ( data%out > 0 ) WRITE( data%out,                                     &
        "( ' ** SUBROUTINE CSGREH: error flag raised during SIF evaluation' )" )
      status = 3
      RETURN

!  non-executable statements

!2000 FORMAT( ' ** SUBROUTINE CSGREH: Increase the size of WK' )
!2020 FORMAT( ' ** SUBROUTINE CSGREH: Increase the size of',                   &
!             ' HE_row_ptr, HE_val_ptr, HE_row or HE_val' )

!  end of subroutine CSGREH

      END SUBROUTINE CSGREH
