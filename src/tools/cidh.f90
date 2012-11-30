! THIS VERSION: CUTEST 1.0 - 24/11/2012 AT 14:40 GMT.

!-*-*-*-*-*-*-*-  C U T E S T    C I D H    S U B R O U T I N E  -*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 77 version originally released in CUTEr, November 1998
!   fortran 2003 version released in CUTEst, 24th November 2012

      SUBROUTINE CIDH( data, status, n, X, iprob, lh1, H )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      TYPE ( CUTEST_data_type ), INTENT( INOUT ) :: data
      INTEGER, INTENT( IN ) :: n, iprob, lh1
      INTEGER, INTENT( OUT ) :: status
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: X
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( lh1, n ) :: H

!  ------------------------------------------------------------
!  compute the Hessian matrix of a specified problem function 
! (iprob = 0 is the objective function, while iprob > 0 is the 
!   iprob-th constraint) of a problem initially written in 
!  Standard Input Format (SIF).

!  H is a two-dimensional array which gives the value of the 
!    Hessian matrix of the proble function evaluated at X. 
!    The i,j-th component of the array will contain the
!    derivative with respect to variables X(i) and X(j)
!  ------------------------------------------------------------

!  local variables

      INTEGER :: i, ig, j, ncalcf, ncalcg, ifstat, igstat, k, nnzh, alloc_status
      REAL ( KIND = wp ) :: ftt
      CHARACTER ( LEN = 80 ) :: bad_alloc = REPEAT( ' ', 80 )
      EXTERNAL :: RANGE 

!  check input parameters

      IF ( lh1 < n ) THEN
        IF ( data%out > 0 ) WRITE( data%out, 2040 ) n
        status = 2 ; RETURN
      END IF
      IF ( iprob < 0 ) THEN
        WRITE( data%out, 2020 ) iprob
        status = 2 ; RETURN
      END IF

!  find group index ig of constraint iprob

      IF ( iprob > 0 ) THEN
        ig = 0
        DO i = 1, data%ng
          IF ( data%KNDOFC( i ) == iprob ) THEN
            ig = i
            EXIT
          END IF
        END DO
        IF ( ig == 0 ) THEN
          WRITE( data%out, 2020 ) iprob
          status = 2 ; RETURN
        END IF
      END IF

!  find which elements are involved in the required problem function.
!  Initialize the list to zero

      DO i = 1, data%nel
        data%ICALCF( i ) = 0
      END DO

!  if group ig is involved, record its elements

      DO ig = 1, data%ng
        IF ( data%KNDOFC( ig ) == iprob ) THEN
          DO j = data%ISTADG( ig ), data%ISTADG( ig + 1 ) - 1
            data%ICALCF( data%IELING( j ) ) = 1
          END DO
        END IF
      END DO

!  only compute the elements which are involved in the required function

      ncalcf = 0
      DO i = 1, data%nel
        IF ( data%ICALCF( i ) == 1 ) THEN
          ncalcf = ncalcf + 1
          data%ICALCF( ncalcf ) = i
        ELSE

!  if this is the first ever evaluation, initialize data%FUVALS

          IF ( data%firstg ) THEN
            data%FUVALS( i ) = 0.0_wp
            DO j = data%INTVAR( i ), data%INTVAR( i + 1 ) - 1
               data%FUVALS( j ) = 0.0_wp
            END DO
            DO j = data%ISTADH( i ), data%ISTADH( i + 1 ) - 1
               data%FUVALS( j ) = 0.0_wp
            END DO
          END IF
        END IF
      END DO

!  evaluate the element function values.

      CALL ELFUN( data%FUVALS, X, data%EPVALU, ncalcf, data%ITYPEE,            &
                  data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,          &
                  data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,           &
                  data%lelvar, data%lntvar, data%lstadh, data%lstep,           &
                  data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,          &
                  1, ifstat )
      IF ( ifstat /= 0 ) GO TO 930

!  evaluate the element function derivatives.

      CALL ELFUN( data%FUVALS, X, data%EPVALU, ncalcf, data%ITYPEE,            &
                  data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,          &
                  data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,           &
                  data%lelvar, data%lntvar, data%lstadh, data%lstep,           &
                  data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,          &
                  3, ifstat )
      IF ( ifstat /= 0 ) GO TO 930

!  compute the list of groups involved in the required problem function.

      ncalcg = 0
      DO ig = 1, data%ng
        IF ( data%KNDOFC( ig ) == iprob ) THEN
          ncalcg = ncalcg + 1
          data%ICALCF( ncalcg ) = ig

!  compute the group argument values ft.

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

!  record the derivatives of trivial groups

          IF ( data%GXEQX( ig ) ) THEN
            data%GVALS( ig, 2 ) = 1.0_wp
            data%GVALS( ig, 3 ) = 0.0_wp
          END IF
        ELSE

!  if this is the first ever evaluation, initialize GVALS

          IF ( data%firstg ) THEN
            data%GVALS( ig, 2 ) = 1.0_wp
            data%GVALS( ig, 3 ) = 0.0_wp
          END IF
        END IF
      END DO

!  evaluate the group derivative values.

      IF ( .NOT. data%altriv ) THEN
        CALL GROUP( data%GVALS, data%ng, data%FT, data%GPVALU, ncalcg,         &
                    data%ITYPEG, data%ISTGP, data%ICALCF, data%ltypeg,         &
                    data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,         &
                    .TRUE., igstat )
        IF ( igstat /= 0 ) GO TO 930
      END IF

!  change the group weightings to isolate the required function

      IF ( data%numcon > 0 ) THEN
        DO ig = 1, data%ng
          i = data%KNDOFC( ig )
          IF ( i == iprob ) THEN
            data%GSCALE_used( ig ) = data%GSCALE( ig )
          ELSE
            data%GSCALE_used( ig ) = 0.0_wp
          END IF
        END DO

!  compute the gradient value

        CALL CUTEST_form_gradients( n, data%ng, data%nel, data%ntotel,         &
               data%nvrels, data%nnza, data%nvargp, data%firstg, data%ICNA,    &
               data%ISTADA, data%IELING, data%ISTADG, data%ISTAEV,             &
               data%IELVAR, data%INTVAR, data%A, data%GVALS( : , 2 ),          &
               data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ),        &
               data%GSCALE_used, data%ESCALE, data%FUVALS( data%lgrjac + 1 ),  &
               data%GXEQX, data%INTREP, data%ISVGRP, data%ISTAGV, data%ITYPEE, &
               data%ISTAJC, data%W_ws, data%W_el, RANGE )
      ELSE
        CALL CUTEST_form_gradients( n, data%ng, data%nel, data%ntotel,         &
               data%nvrels, data%nnza, data%nvargp, data%firstg, data%ICNA,    &
               data%ISTADA, data%IELING, data%ISTADG, data%ISTAEV,             &
               data%IELVAR, data%INTVAR, data%A, data%GVALS( : , 2 ),          &
               data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ),        &
               data%GSCALE, data%ESCALE, data%FUVALS( data%lgrjac + 1 ),       &
               data%GXEQX, data%INTREP, data%ISVGRP, data%ISTAGV, data%ITYPEE, &
               data%ISTAJC, data%W_ws, data%W_el, RANGE )
      END IF

!      CALL ELGRD( n, data%ng, data%firstg, data%ICNA( 1 ), data%licna, &
!                      data%ISTADA( 1 ), data%lstada, data%IELING( 1 ), &
!                      data%leling, data%ISTADG( 1 ), data%lstadg, &
!                      data%ITYPEE( 1 ), data%lintre, &
!                      data%ISTAEV( 1 ), data%lstaev, data%IELVAR( 1 ), &
!                      data%lelvar, data%INTVAR( 1 ), data%lntvar, &
!                      data%IWORK( data%lsvgrp + 1 ), &
!                      data%lnvgrp, data%IWORK( data%lstajc + 1 ), data%lnstjc, &
!                      data%IWORK( data%lstagv + 1 ), data%lnstgv, &
!                      data%A( 1 ), data%la, &
!                      data%GVALS(  : , 2 ), data%lgvals, &
!                      data%FUVALS, data%lnguvl, &
!                      data%FUVALS( data%lggfx + 1 ), &
!                      data%WRK( 1 ), data%ng, &
!                      data%ESCALE( 1 ), data%lescal, &
!                      data%FUVALS( data%lgrjac + 1 ), &
!                      data%lngrjc, data%WRK( 1 ), data%WRK( n + 1 ), &
!                      data%maxsel, &
!                      data%GXEQX( 1 ), data%lgxeqx, &
!                      data%INTREP( 1 ), data%lintre, RANGE )
!      ELSE
!
!  Compute the gradient value
!
!      CALL ELGRD( n, data%ng, data%firstg, data%ICNA( 1 ), data%licna, &
!                      data%ISTADA( 1 ), data%lstada, data%IELING( 1 ), &
!                      data%leling, data%ISTADG( 1 ), data%lstadg, &
!                      data%ITYPEE( 1 ), data%lintre, &
!                      data%ISTAEV( 1 ), data%lstaev, data%IELVAR( 1 ), &
!                      data%lelvar, data%INTVAR( 1 ), data%lntvar, &
!                      data%IWORK( data%lsvgrp + 1 ), &
!                      data%lnvgrp, data%IWORK( data%lstajc + 1 ), &
!                      data%lnstjc, &
!                      data%IWORK( data%lstagv + 1 ), data%lnstgv, &
!                      data%A( 1 ), data%la, &
!                      data%GVALS( : , 2 ), data%lgvals, &
!                      data%FUVALS, data%lnguvl, &
!                      data%FUVALS( data%lggfx + 1 ), &
!                      data%GSCALE( 1 ), data%lgscal, &
!                      data%ESCALE( 1 ), data%lescal, &
!                      data%FUVALS( data%lgrjac + 1 ), &
!                      data%lngrjc, data%WRK( 1 ), data%WRK( n + 1 ), &
!                      data%maxsel, &
!                      data%GXEQX( 1 ), data%lgxeqx, &
!                      data%INTREP( 1 ), data%lintre, RANGE )

      data%firstg = .FALSE.

!  assemble the Hessian; use every variable

!      DO i = 1, n
!        data%IVAR( i ) = i
!      END DO

      IF ( data%numcon > 0 ) THEN
        CALL CUTEST_assemble_hessian(                                          &
             n, data%ng, data%nel, data%ntotel, data%nvrels, data%nnza,        &
             data%maxsel, data%nvargp, data%ISTADH,                            &
             data%ICNA, data%ISTADA, data%INTVAR, data%IELVAR, data%IELING,    &
             data%ISTADG, data%ISTAEV, data%ISTAGV, data%ISVGRP, data%A,       &
             data%FUVALS, data%lnguvl, data%FUVALS, data%lnhuvl,               &
             data%GVALS( : , 2 ), data%GVALS( :  , 3 ), data%GSCALE_used,      &
             data%ESCALE, data%GXEQX, data%ITYPEE, data%INTREP, RANGE,         &
             0, data%out, data%out, data%io_buffer, .FALSE., .FALSE.,          &
             n, status, alloc_status, bad_alloc,                               &
             data%array_status, data%lh_row, data%lh_col, data%lh_val,         &
             data%H_row, data%H_col, data%H_val,                               &
             data%LINK_col, data%POS_in_H, data%llink, data%lpos,              &
             data%W_ws, data%W_el, data%W_in, data%H_el, data%H_in,            &
             nnzh = nnzh )
      ELSE
        CALL CUTEST_assemble_hessian(                                          &
             n, data%ng, data%nel, data%ntotel, data%nvrels, data%nnza,        &
             data%maxsel, data%nvargp, data%ISTADH,                            &
             data%ICNA, data%ISTADA, data%INTVAR, data%IELVAR, data%IELING,    &
             data%ISTADG, data%ISTAEV, data%ISTAGV, data%ISVGRP, data%A,       &
             data%FUVALS, data%lnguvl, data%FUVALS, data%lnhuvl,               &
             data%GVALS( : , 2 ), data%GVALS( :  , 3 ), data%GSCALE,           &
             data%ESCALE, data%GXEQX, data%ITYPEE, data%INTREP, RANGE,         &
             0, data%out, data%out, data%io_buffer, .FALSE., .FALSE.,          &
             n, status, alloc_status, bad_alloc,                               &
             data%array_status, data%lh_row, data%lh_col, data%lh_val,         &
             data%H_row, data%H_col, data%H_val,                               &
             data%LINK_col, data%POS_in_H, data%llink, data%lpos,              &
             data%W_ws, data%W_el, data%W_in, data%H_el, data%H_in,            &
             nnzh = nnzh )
      END IF

!  check for errors in the assembly

      IF ( status > 0 ) RETURN

!!$      IF ( data%numcon > 0 ) THEN
!!$      CALL ASMBL( n, data%ng, data%maxsel, n, lh, lih, nnzh, &
!!$                   n, data%IVAR( 1), data%ISTADH( 1 ), data%lstadh, &
!!$                   data%ICNA( 1 ), data%licna, &
!!$                   data%ISTADA( 1 ), data%lstada, &
!!$                   data%INTVAR( 1 ), data%lntvar, &
!!$                   data%IELVAR( 1 ), data%lelvar, &
!!$                   data%IELING( 1 ), data%leling, &
!!$                   data%ISTADG( 1 ), data%lstadg, &
!!$                   data%ISTAEV( 1 ), data%lstaev, &
!!$                   data%IWORK( data%lstagv + 1 ), data%lnstgv, &
!!$                   data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, &
!!$                   data%IWORK( data%lsend + lirnh + 1 ), &
!!$                   data%IWORK( data%lsend + ljcnh + 1 ), &
!!$                   data%IWORK( data%lsend + lnxtrw + 1 ), linxtr, &
!!$                   data%IWORK( data%lsend + liwkh + 1 ), n, &
!!$                   data%A( 1 ), data%la, data%FUVALS, data%lnguvl, &
!!$                   data%FUVALS, data%lnhuvl, &
!!$                   data%GVALS( : , 2 ), data%GVALS( : , 3 ), &
!!$                   data%WRK( 1 ), data%ESCALE( 1 ), data%lescal, &
!!$                   data%WRK( data%ng + 1 ), data%WRK( lwkh + data%ng + 1 ), &
!!$                   data%lwk2 - lwkh, data%GXEQX( 1 ), &
!!$                   data%lgxeqx, data%INTREP( 1 ), data%lintre, &
!!$                   data%ITYPEE( 1 ), data%lintre, &
!!$                   RANGE, 1, data%out, .FALSE., i, inform, .FALSE., &
!!$                   .FALSE. )
!!$      ELSE
!!$      CALL ASMBL( n, data%ng, data%maxsel, n, lh, lih, nnzh, &
!!$                   n, data%IVAR( 1), data%ISTADH( 1 ), data%lstadh, &
!!$                   data%ICNA( 1 ), data%licna, &
!!$                   data%ISTADA( 1 ), data%lstada, &
!!$                   data%INTVAR( 1 ), data%lntvar, &
!!$                   data%IELVAR( 1 ), data%lelvar, &
!!$                   data%IELING( 1 ), data%leling, &
!!$                   data%ISTADG( 1 ), data%lstadg, &
!!$                   data%ISTAEV( 1 ), data%lstaev, &
!!$                   data%IWORK( data%lstagv + 1 ), &
!!$                   data%lnstgv, data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, &
!!$                   data%IWORK( data%lsend + lirnh + 1 ), &
!!$                   data%IWORK( data%lsend + ljcnh + 1 ), &
!!$                   data%IWORK( data%lsend + lnxtrw + 1 ), linxtr, &
!!$                   data%IWORK( data%lsend + liwkh + 1 ), n, &
!!$                   data%A( 1 ), data%la, &
!!$                   data%FUVALS, data%lnguvl, data%FUVALS, data%lnhuvl, &
!!$                   data%GVALS( : , 2 ), data%GVALS( : , 3 ), &
!!$                   data%GSCALE( 1 ), data%ESCALE( 1 ), data%lescal, &
!!$                   data%WRK( 1 ), data%WRK( lwkh + 1 ), &
!!$                   data%lwk2 - lwkh, data%GXEQX( 1 ), &
!!$                   data%lgxeqx, data%INTREP( 1 ), data%lintre, &
!!$                   data%ITYPEE( 1 ), data%lintre, &
!!$                   RANGE, 1, data%out, .FALSE., i, inform, .FALSE., &
!!$                   .FALSE. )
!!$      END IF

!  initialize the dense matrix

      H( : n, : n ) = 0.0_wp

!  transfer the matrix from co-ordinate to dense storage and symmetrize the 
!  martix

      DO k = 1, nnzh
        i = data%H_row( k ) ; j = data%H_col( k )
        H( i, j ) = data%H_val( k ) ; H( j, i ) = data%H_val( k )
      END DO

!!$      IF ( data%numcon > 0 ) THEN
!!$        DO k = 1, nnzh
!!$          i = data%IWORK( data%lsend + lirnh + k )
!!$          j = data%IWORK( data%lsend + ljcnh + k )
!!$          H( i, j ) = data%WRK( data%ng + k )
!!$          H( j, i ) = data%WRK( data%ng + k )
!!$        END DO
!!$      ELSE
!!$        DO k = 1, nnzh
!!$          i = data%IWORK( data%lsend + lirnh + k )
!!$          j = data%IWORK( data%lsend + ljcnh + k )
!!$          H( i, j ) = data%WRK( k )
!!$          H( j, i ) = data%WRK( k )
!!$        END DO
!!$      END IF

!  update the counters for the report tool

      IF ( iprob == 0 ) THEN
        data%nc2oh = data%nc2oh + 1
      ELSE 
        data%nc2ch = data%nc2ch + 1
      ENDIF 
      status = 0
      RETURN

!  unsuccessful returns

  930 CONTINUE
      IF ( data%out > 0 ) WRITE( data%out,                                     &
        "( ' ** SUBROUTINE CIDH: error flag raised during SIF evaluation' )" )
      status = 3
      RETURN

!  non-executable statements

 2020 FORMAT( ' ** SUBROUTINE CIDH: invalid problem index iprob = ', I0 )
 2040 FORMAT( ' ** SUBROUTINE CIDH: Increase the leading dimension of H to ',  &
              I0 )

!  end of subroutine CIDH

      END SUBROUTINE CIDH
