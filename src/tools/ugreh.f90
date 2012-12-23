! THIS VERSION: CUTEST 1.0 - 27/11/2012 AT 14:00 GMT.

!-*-*-*-*-*-*-  C U T E S T    U G R E H    S U B R O U T I N E  -*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 77 version originally released in CUTEr, November 1994
!   fortran 2003 version released in CUTEst, 27th November 2012

      SUBROUTINE CUTEST_ugreh( data, work, status, n, X, G, ne, lhe_ptr, HE_row_ptr, &
                        HE_val_ptr, lhe_row, HE_row, lhe_val, HE_val, byrows )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      TYPE ( CUTEST_data_type ), INTENT( IN ) :: data
      TYPE ( CUTEST_work_type ), INTENT( INOUT ) :: work
      INTEGER, INTENT( IN ) :: n, lhe_ptr, lhe_row, lhe_val
      INTEGER, INTENT( OUT ) :: ne, status
      LOGICAL, INTENT( IN ) :: byrows
      INTEGER, INTENT( OUT ), DIMENSION( lhe_ptr ) :: HE_row_ptr, HE_val_ptr
      INTEGER, INTENT( OUT ), DIMENSION( lhe_row ) :: HE_row
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: X
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: G
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( lhe_val ) :: HE_val

!  ----------------------------------------------------------------------------
!  compute the gradient and Hessian matrix of a group partially separable 
!  function initially written in Standard Input Format (SIF)

!  the matrix is represented in "finite element format", i.e., 

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

!  local variables

      INTEGER :: i, ig, j, ifstat, igstat, lhe_row_int, lhe_val_int
      INTEGER :: alloc_status
      REAL ( KIND = wp ) :: ftt
      CHARACTER ( LEN = 80 ) :: bad_alloc = REPEAT( ' ', 80 )
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

!  evaluate the element function gradients and Hessians

      CALL ELFUN( work%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,          &
                  data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,          &
                  data%ISTEP, work%ICALCF, data%ltypee, data%lstaev,           &
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
          ftt = ftt + data%ESCALE( j ) * work%FUVALS( data%IELING( j ) )
        END DO
        work%FT( ig ) = ftt

!  record the derivatives of trivial groups

        IF ( data%GXEQX( ig ) ) THEN
          work%GVALS( ig, 2 ) = 1.0_wp
          work%GVALS( ig, 3 ) = 0.0_wp
        END IF
      END DO

!  evaluate the group derivative values

      IF ( .NOT. data%altriv ) THEN
        CALL GROUP( work%GVALS, data%ng, work%FT, data%GPVALU, data%ng,        &
                    data%ITYPEG, data%ISTGP, work%ICALCF, data%ltypeg,         &
                    data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,         &
                    .TRUE., igstat )
        IF ( igstat /= 0 ) GO TO 930
      END IF

!  compute the gradient value

      CALL CUTEST_form_gradients( n, data%ng, data%nel, data%ntotel,           &
             data%nvrels, data%nnza, data%nvargp, work%firstg, data%ICNA,      &
             data%ISTADA, data%IELING, data%ISTADG, data%ISTAEV,               &
             data%IELVAR, data%INTVAR, data%A, work%GVALS( : , 2 ),            &
             work%FUVALS, data%lnguvl, work%FUVALS( data%lggfx + 1 ),          &
             data%GSCALE, data%ESCALE, work%FUVALS( data%lgrjac + 1 ),         &
             data%GXEQX, data%INTREP, data%ISVGRP, data%ISTAGV, data%ITYPEE,   &
             work%ISTAJC, work%W_ws, work%W_el, RANGE )
      work%firstg = .FALSE.

!  store the gradient value

      DO i = 1, n
        G( i ) = work%FUVALS( data%lggfx + i )
      END DO

!  assemble the Hessian

      lhe_row_int = lhe_row ; lhe_val_int = lhe_val
      CALL CUTEST_assemble_element_hessian(                                    &
                        data%ng, data%nel,data% ntotel, data%nvrels,           &
                        data%nnza, data%maxsel, data%nvargp,                   &
                        data%lnguvl, data%lnhuvl, data%ISTADH, data%ICNA,      &
                        data%ISTADA, data%INTVAR, data%IELVAR,                 &
                        data%IELING, data%ISTADG, data%ISTAEV,                 &
                        data%ISTAGV, data%ISVGRP, data%ITYPEE,                 &
                        data%A, work%FUVALS, work%FUVALS,                      &
                        work%GVALS( : , 2 ), work%GVALS( : , 3 ),              &
                        data%GSCALE, data%ESCALE, data%GXEQX, data%INTREP,     &
                        work%ISWKSP, work%W_ws, work%W_el, work%W_in,          &
                        work%H_el, work%H_in, RANGE, ne, lhe_row_int,          &
                        lhe_val_int, work%H_row, HE_row_ptr, work%H_val,       &
                        HE_val_ptr, byrows, 0, data%out, data%out,             &
                        data%io_buffer, alloc_status, bad_alloc, status )

!  check for errors in the assembly

      IF ( status > 0 ) RETURN

!  check that HE_row and HE_val are large enough

      IF ( lhe_row < HE_row_ptr( ne + 1 ) - 1 ) THEN
        IF ( data%out > 0 ) WRITE( data%out, "( ' ** SUBROUTINE UGREH: ',      &
       &  'Increase the dimension of HE_row to ',  I0 )" )                     &
             HE_row_ptr( ne + 1 ) - 1 
        status = 2 ; RETURN
      END IF

      IF ( lhe_val < HE_val_ptr( ne + 1 ) - 1 ) THEN
        IF ( data%out > 0 ) WRITE( data%out, "( ' ** SUBROUTINE UGREH: ',      &
       &  'Increase the dimension of HE_val to ',  I0 )" )                     &
             HE_val_ptr( ne + 1 ) - 1 
        status = 2 ; RETURN
      END IF

!  record the element Hessian

      HE_row( : HE_row_ptr( ne + 1 ) - 1 )                                     &
         = work%H_row( : HE_row_ptr( ne + 1 ) - 1 )
      HE_val( : HE_val_ptr( ne + 1 ) - 1 )                                     &
         = work%H_val( : HE_val_ptr( ne + 1 ) - 1 )

!  update the counters for the report tool

      work%nc2og = work%nc2og + 1
      work%nc2oh = work%nc2oh + 1
      status = 0
      RETURN

!  unsuccessful returns

  930 CONTINUE
      IF ( data%out > 0 ) WRITE( data%out,                                     &
        "( ' ** SUBROUTINE UGREH: error flag raised during SIF evaluation' )" )
      status = 3
      RETURN

!  end of subroutine CUTEST_ugreh

      END SUBROUTINE CUTEST_ugreh
