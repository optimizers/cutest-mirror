! THIS VERSION: CUTEST 1.0 - 23/11/2012 AT 15:00 GMT.

!-*-*-*-*-*-*-*-  C U T E S T    U S H    S U B R O U T I N E  -*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 77 version originally released in CUTE, July 1991
!   fortran 2003 version released in CUTEst, 23rd November 2012

      SUBROUTINE USH( data, status, n, X, nnzh, lh, H_val, H_row, H_col )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      TYPE ( CUTEST_data_type ), INTENT( INOUT ) :: data
      INTEGER, INTENT( IN ) :: n, lh
      INTEGER, INTENT( OUT ) :: nnzh, status
      INTEGER, INTENT( OUT ), DIMENSION( lh ) :: H_row, H_col
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: X
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( lh ) :: H_val

!  -------------------------------------------------------------------
!  compute the Hessian matrix of a group partially separable function.
!  The upper triangle of the Hessian is stored in coordinate form,
!  i.e., the entry H_val(i) has row index H_row(i) and column index 
!  H_col(i) for i = 1, ...., nnzh
!  ------------------------------------------------------------------

!  Based on the minimization subroutine data%laNCELOT/SBMIN
!  by Conn, Gould and Toint.

!  Nick Gould, for CGT productions.
!  July 1991.

!  Local variables

      INTEGER :: i, ig, j, ifstat, igstat, alloc_status
      REAL ( KIND = wp ) :: ftt
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

!  evaluate the element function values

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

!  include the contributions from the nonlinear elements

        DO j = data%ISTADG( ig ), data%ISTADG( ig + 1 ) - 1
          ftt = ftt + data%ESCALE( j ) * data%FUVALS( data%IELING( j ) )
        END DO
        data%FT( ig ) = ftt

!  record the derivatives of trivial groups

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

!  compute the gradient value

!      CALL ELGRD( n, data%ng, data%firstg, data%ICNA( 1 ), data%licna, &
!                   data%ISTADA( 1 ), data%lstada, data%IELING( 1 ), &
!                   data%leling, data%ISTADG( 1 ), data%lstadg, &
!                   data%ITYPEE( 1 ), data%lintre, &
!                   data%ISTAEV( 1 ), data%lstaev, &
!                   data%IELVAR( 1 ), data%lelvar, &
!                   data%INTVAR( 1 ), data%lntvar, &
!                   data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, &
!                   data%IWORK( data%lstajc + 1 ), data%lnstjc, &
!                   data%IWORK( data%lstagv + 1 ), data%lnstgv, &
!                   data%A( 1 ), data%la, &
!                   data%GVALS( : , 2 ), data%lgvals, &
!                   data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ), &
!                   data%GSCALE( 1 ), data%lgscal, &
!                   data%ESCALE( 1 ), data%lescal, &
!                   data%FUVALS( data%lgrjac + 1 ), &
!                   data%lngrjc, data%WRK( 1 ), data%WRK( n + 1 ), &
!                   data%maxsel, &
!                   data%GXEQX( 1 ), data%lgxeqx, &
!                   data%INTREP( 1 ), data%lintre, RANGE )

      CALL CUTEST_form_gradients( n, data%ng, data%nel, data%ntotel,           &
             data%nvrels, data%nnza, data%nvargp, data%firstg, data%ICNA,      &
             data%ISTADA, data%IELING, data%ISTADG, data%ISTAEV,               &
             data%IELVAR, data%INTVAR, data%A, data%GVALS( : , 2 ),            &
             data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ),          &
             data%GSCALE, data%ESCALE, data%FUVALS( data%lgrjac + 1 ),         &
             data%GXEQX, data%INTREP, data%ISVGRP, data%ISTAGV, data%ITYPEE,   &
             data%ISTAJC, data%W_ws, data%W_el, RANGE )
      data%firstg = .FALSE.

!  assemble the Hessian; use every variable

      DO i = 1, n
        data%IVAR( i ) = i
      END DO

      CALL CUTEST_assemble_hessian(                                            &
             n, data%ng, data%nel, data%ntotel, data%nvrels, data%nnza,        &
             data%maxsel, data%nvargp, n, data%IVAR, data%ISTADH,              &
             data%ICNA, data%ISTADA, data%INTVAR, data%IELVAR, data%IELING,    &
             data%ISTADG, data%ISTAEV, data%ISTAGV, data%ISVGRP, data%A,       &
             data%FUVALS, data%lnguvl, data%FUVALS, data%lnhuvl,               &
             data%GVALS( : , 2 ), data%GVALS( :  , 3 ), data%GSCALE,           &
             data%ESCALE, data%GXEQX, data%ITYPEE, data%INTREP, RANGE,         &
             0, data%out, data%out, data%io_buffer, .FALSE.,                   &
             .FALSE., .FALSE., 0, status, alloc_status, bad_alloc,             &
             data%assemble_data, data%lirnh, data%ljcnh, data%lh,              &
             data%H_row, data%H_col, data%H_val,                               &
             data%LINK_col, data%POS_in_H, data%llink, data%lpos,              &
             data%IW_asmbl, data%W_ws, data%W_el, data%W_in, data%H_el,        &
             data%H_in, data%skipg, nnzh = nnzh )

      H_row( : nnzh ) = data%H_row( : nnzh )
      H_col( : nnzh ) = data%H_col( : nnzh )
      H_val( : nnzh ) = data%H_val( : nnzh )

!      CALL ASMBL( n, data%ng, data%maxsel, n, lh, lih, nnzh, &
!                   n, data%IVAR( 1), data%ISTADH( 1 ), data%lstadh, &
!                   data%ICNA( 1 ), data%licna, &
!                   data%ISTADA( 1 ), data%lstada, &
!                   data%INTVAR( 1 ), data%lntvar, &
!                   data%IELVAR( 1 ), data%lelvar, &
!                   data%IELING( 1 ), data%leling, &
!                   data%ISTADG( 1 ), data%lstadg, &
!                   data%ISTAEV( 1 ), data%lstaev, &
!                   data%IWORK( data%lstagv + 1 ), data%lnstgv, &
!                   data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, &
!                   IRNH, ICNH, data%IWORK( data%lsend + lnxtrw + 1 ), linxtr, &
!                   data%IWORK( data%lsend + liwkh + 1 ), n, &
!                   data%A( 1 ), data%la, &
!                   data%FUVALS, data%lnguvl, data%FUVALS, data%lnhuvl, &
!                   data%GVALS( : , 2 ), data%GVALS( : , 3 ), &
!                   data%GSCALE( 1 ), data%ESCALE( 1 ), data%lescal, &
!                   H, data%WRK( 1 ), data%lwk2, data%GXEQX( 1 ), &
!                   data%lgxeqx, data%INTREP( 1 ), data%lintre, &
!                   data%ITYPEE( 1 ), data%lintre, &
!                   RANGE, 1, data%out, .FALSE., i, inform, .FALSE., &
!                  .TRUE. )

!  check for errors in the assembly

      IF ( status > 0 ) RETURN

!  update the counters for the report tool

      data%nc2oh = data%nc2oh + 1
      status = 0
      RETURN

!  unsuccessful returns

  930 CONTINUE
      IF ( data%out > 0 ) WRITE( data%out,                                     &
        "( ' ** SUBROUTINE USH: error flag raised during SIF evaluation' )" )
      status = 3
      RETURN

!  end of subroutine USH

      END SUBROUTINE USH
