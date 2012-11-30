! THIS VERSION: CUTEST 1.0 - 24/11/2012 AT 14:30 GMT.

!-*-*-*-*-*-*-  C U T E S T    C G R D H    S U B R O U T I N E  -*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 77 version originally released in CUTE, November 1991
!   fortran 2003 version released in CUTEst, 24th November 2012

      SUBROUTINE CGRDH( data, status, n, m, X, Y, grlagf, G, jtrans,           &
                        lcjac1, lcjac2, CJAC, lh1, H )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      TYPE ( CUTEST_data_type ), INTENT( INOUT ) :: data
      INTEGER, INTENT( IN ) :: n, m, lcjac1, lcjac2, lh1
      INTEGER, INTENT( OUT ) :: status
      LOGICAL, INTENT( IN ) :: grlagf, jtrans
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: X
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( m ) :: Y
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: G
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( lcjac1, lcjac2 ) :: CJAC
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( lh1, n ) :: H

!  ---------------------------------------------------------------
!  compute both the gradients of the objective, or Lagrangian, and
!  general constraint functions and the Hessian matrix of the
!  Lagrangian function of a problem initially written in Standard
!  Input Format (SIF)

!  G	 is an array which gives the value of the gradient of
!	 the objective function evaluated at X (grlagf = .FALSE.)
!        of of the Lagrangian function evaluated at X and Y
!        (GRLAGF = .TRUE.)

!  CJAC	 is a two-dimensional array of dimension ( lcjac1, lcjac2 )
!	 which gives the value of the Jacobian matrix of the
!	 constraint functions, or its transpose, evaluated at X.
!	 If jtrans is .TRUE., the i,j-th component of the array
!        will contain the i-th derivative of the j-th constraint
!        function. Otherwise, if jtrans is .FALSE., the i,j-th
!        component of the array will contain the j-th derivative
!        of the i-th constraint function

!  H     is a two-dimensional array which gives the value of the
!        Hessian matrix of the Lagrangian function evaluated at
!        X and Y. The i,j-th component of the array will contain 
!        the derivative with respect to variables X(i) and X(j)
!  ---------------------------------------------------------------

!  local variables

      INTEGER :: i, j, iel, k, ig, ii, ig1, l, jj, ll, nnzh, iendgv, icon
      INTEGER :: nin, nvarel, nelow, nelup, istrgv, ifstat, igstat, alloc_status
      REAL ( KIND = wp ) :: ftt, gi, gii, scalee
      CHARACTER ( LEN = 80 ) :: bad_alloc = REPEAT( ' ', 80 )
      LOGICAL :: nontrv
      EXTERNAL :: RANGE 

!  check input parameters

      IF ( data%numcon > 0 ) THEN
        IF ( jtrans ) THEN
          IF ( lcjac1 < n .OR. lcjac2 < m ) THEN
            IF ( lcjac1 < n .AND. data%out > 0 ) WRITE( data%out, 2020 ) n
            IF ( lcjac2 < m .AND. data%out > 0 ) WRITE( data%out, 2030 ) m
            status = 2 ; RETURN
          END IF
        ELSE
          IF ( lcjac1 < m .OR. lcjac2 < n ) THEN
            IF ( lcjac1 < m .AND. data%out > 0 ) WRITE( data%out, 2020 ) m
            IF ( lcjac2 < n .AND. data%out > 0 ) WRITE( data%out, 2030 ) n
            status = 2 ; RETURN
          END IF
        END IF
      END IF
      IF ( lh1 < n ) THEN
        IF ( data%out > 0 ) WRITE( data%out, 2040 ) n
        status = 2 ; RETURN
      END IF

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

!  for unconstrained problems, skip construction of data%W_ws( ig ) and skip 
!  specialized construction of gradient and Jacobian.  Call ELGRD instead

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

!  compute the gradient values. Initialize the gradient and Jacobian (or its 
!  transpose) as zero

        G( : n ) = 0.0_wp
        IF ( jtrans ) THEN
          CJAC( : n, : m ) = 0.0_wp
        ELSE
          CJAC( : m, : n ) = 0.0_wp
        END IF

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
          gii = data%GSCALE_used( ig )
          IF  ( nontrv ) THEN
            gi = gi * data%GVALS( ig, 2 )
            gii = gii * data%GVALS( ig, 2 )
          END IF

!  this is the first gradient evaluation or the group has nonlinear elements

          IF ( data%firstg .OR. nelow <= nelup ) THEN
            data%W_ws( data%ISVGRP( istrgv : iendgv ) ) = 0.0_wp

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
                             data%W_el, nvarel, nin, data%ITYPEE( iel ),       &
                             nin, nvarel )
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

!  Allocate a gradient

!DIR$ IVDEP
             DO i = istrgv, iendgv
               ll = data%ISVGRP( i )

!  the group belongs to the objective function

               IF ( icon == 0 ) THEN
                 G( ll ) = G( ll ) + gi * data%W_ws( ll )

!  the group defines a constraint

               ELSE
                 IF ( JTRANS ) THEN
                   CJAC( ll, icon ) = gi * data%W_ws( ll )
                 ELSE
                   CJAC( icon, ll ) = gi * data%W_ws( ll )
                 END IF
                 IF ( grlagf ) G( ll ) = G( ll ) + gii * data%W_ws( ll )
               END IF

!  if the group is non-trivial, also store the nonzero entries of the
!  gradient of the function in GRJAC

               IF ( nontrv ) THEN
                 jj = data%ISTAJC( ll )
                 data%FUVALS( data%lgrjac + jj ) = data%W_ws( ll )

!  increment the address for the next nonzero in the column of the Jacobian 
!  for variable ll

                 data%ISTAJC( ll ) = jj + 1
               END IF
             END DO

!  This is not the first gradient evaluation and there is only a linear element

          ELSE

!  allocate a gradient

!DIR$ IVDEP
            DO k = data%ISTADA( ig ), data%ISTADA( ig1 ) - 1
              ll = data%ICNA( k )

!  the group belongs to the objective function

              IF ( icon == 0 ) THEN
                G( ll ) = G( ll ) + gi * data%A( k )

!  the group defines a constraint

              ELSE
                IF ( JTRANS ) THEN
                  CJAC( ll, icon ) = gi * data%A( k )
                ELSE
                  CJAC( icon, ll ) = gi * data%A( k )
                END IF
                IF ( grlagf ) G( ll ) = G( ll ) + gii * data%A( k )
              END IF
            END DO

!  the group is non-trivial; increment the starting addresses for the groups 
!  used by each variable in the (unchanged) linear element to avoid resetting
!  the nonzeros in the Jacobian

            IF ( nontrv ) THEN
!DIR$ IVDEP
              DO i = istrgv, iendgv
                ll = data%ISVGRP( i )
                data%ISTAJC( ll ) = data%ISTAJC( ll ) + 1
              END DO
            END IF
          END IF
        END DO

!  reset the starting addresses for the lists of groups using each variable to 
!  their values on entry

        DO i = n, 2, - 1
           data%ISTAJC( i ) = data%ISTAJC( i - 1 )
        END DO
        data%ISTAJC( 1 ) = 1
      ELSE

!  compute the gradient value

!        CALL ELGRD( n, data%ng, data%firstg, data%ICNA( 1 ), data%licna, &
!                      data%ISTADA( 1 ), data%lstada, data%IELING( 1 ), &
!                      data%leling, data%ISTADG( 1 ), data%lstadg, &
!                      data%ITYPEE( 1 ), data%lintre, &
!                      data%ISTAEV( 1 ), data%lstaev, data%IELVAR( 1 ), &
!                      data%lelvar, data%INTVAR( 1 ), data%lntvar, &
!                      data%ISVGRP( 1 ), data%lnvgrp, &
!                      data%ISTAJC( 1 ), data%lnstjc, &
!                      data%ISTAGV( 1 ), data%lnstgv, &
!                      data%A( 1 ), data%la, &
!                      data%GVALS( : , 2 ), data%lgvals, &
!                      data%FUVALS, data%lnguvl, &
!                      data%FUVALS( data%lggfx + 1 ), &
!                      data%GSCALE( 1 ), data%lgscal, &
!                      data%ESCALE( 1 ), data%lescal, &
!                      data%FUVALS( data%lgrjac + 1 ), &
!                      data%lngrjc, data%W_ws( 1 ), data%W_el( 1 ), &
!                      data%maxsel, &
!                      data%GXEQX( 1 ), data%lgxeqx, &
!                      data%INTREP( 1 ), data%lintre, RANGE )

        CALL CUTEST_form_gradients( n, data%ng, data%nel, data%ntotel,         &
               data%nvrels, data%nnza, data%nvargp, data%firstg, data%ICNA,    &
               data%ISTADA, data%IELING, data%ISTADG, data%ISTAEV,             &
               data%IELVAR, data%INTVAR, data%A, data%GVALS( : , 2 ),          &
               data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ),        &
               data%GSCALE, data%ESCALE, data%FUVALS( data%lgrjac + 1 ),       &
               data%GXEQX, data%INTREP, data%ISVGRP, data%ISTAGV, data%ITYPEE, &
               data%ISTAJC, data%W_ws, data%W_el, RANGE )

!  store the gradient value

        DO i = 1, n
          G( i ) = data%FUVALS( data%lggfx + i )
        END DO
      END IF
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
             0, data%out, data%out, data%io_buffer, .TRUE., .FALSE.,           &
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
             0, data%out, data%out, data%io_buffer, .TRUE., .FALSE.,           &
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
!!$                   data%ISTAGV( 1 ), data%lnstgv, &
!!$                   data%ISVGRP( 1 ), data%lnvgrp, &
!!$                   data%IWORK( data%lsend + lirnh + 1 ), &
!!$                   data%IWORK( data%lsend + ljcnh + 1 ), &
!!$                   data%IWORK( data%lsend + lnxtrw + 1 ), linxtr, &
!!$                   data%IWORK( data%lsend + liwkh + 1 ), n, &
!!$                   data%A( 1 ), data%la, &
!!$                   data%FUVALS, data%lnguvl, data%FUVALS, data%lnhuvl, &
!!$                   data%GVALS( : , 2 ), data%GVALS( : , 3 ), &
!!$                   data%W_ws( 1 ), data%ESCALE( 1 ), data%lescal, &
!!$                   data%W_ws( data%ng + 1 ), data%W_ws( lwkh + data%ng + 1 ), &
!!$                   data%lwk2 - lwkh, data%GXEQX( 1 ), &
!!$                   data%lgxeqx, data%INTREP( 1 ), data%lintre, &
!!$                   data%ITYPEE( 1 ), data%lintre, &
!!$                   RANGE, 1, data%out, .FALSE., i, inform, &
!!$                   .FALSE., .FALSE. )
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
!!$                   data%ISTAGV( 1 ), data%lnstgv, &
!!$                   data%ISVGRP( 1 ), data%lnvgrp, &
!!$                   data%IWORK( data%lsend + lirnh + 1 ), &
!!$                   data%IWORK( data%lsend + ljcnh + 1 ), &
!!$                   data%IWORK( data%lsend + lnxtrw + 1 ), linxtr, &
!!$                   data%IWORK( data%lsend + liwkh + 1 ), n, &
!!$                   data%A( 1 ), data%la, &
!!$                   data%FUVALS, data%lnguvl, data%FUVALS, data%lnhuvl, &
!!$                   data%GVALS( : , 2 ), data%GVALS( : , 3 ), &
!!$                   data%GSCALE( 1 ), data%ESCALE( 1 ), data%lescal, &
!!$                   data%W_ws( 1 ), data%W_ws( lwkh + 1 ), &
!!$                   data%lwk2 - lwkh, data%GXEQX( 1 ), &
!!$                   data%lgxeqx, data%INTREP( 1 ), data%lintre, &
!!$                   data%ITYPEE( 1 ), data%lintre, &
!!$                   RANGE, 1, data%out, .FALSE., i, inform, &
!!$                   .FALSE., .FALSE. )
!!$      END IF

!  initialize the dense Hessian matrix

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
!!$          H( i, j ) = data%W_ws( data%ng + k )
!!$          H( j, i ) = data%W_ws( data%ng + k )
!!$        END DO
!!$      ELSE
!!$        DO k = 1, nnzh
!!$          i = data%IWORK( data%lsend + lirnh + k )
!!$          j = data%IWORK( data%lsend + ljcnh + k )
!!$          H( i, j ) = data%W_ws( k )
!!$          H( j, i ) = data%W_ws( k )
!!$        END DO
!!$      END IF

!  update the counters for the report tool

      data%nc2og = data%nc2og + 1
      data%nc2cg = data%nc2cg + data%pnc
      data%nc2oh = data%nc2oh + 1
      data%nc2ch = data%nc2ch + data%pnc
      status = 0
      RETURN

!  unsuccessful returns

  930 CONTINUE
      IF ( data%out > 0 ) WRITE( data%out,                                     &
        "( ' ** SUBROUTINE CGRDH: error flag raised during SIF evaluation' )" )
      status = 3
      RETURN

!  non-executable statements

 2020 FORMAT( ' ** SUBROUTINE CGRDH: Increase the leading dimension of CJAC',  &
              ' to ', I0 )
 2030 FORMAT( ' ** SUBROUTINE CGRDH: Increase the second dimension of CJAC',   &
              ' to ', I0 )
 2040 FORMAT( ' ** SUBROUTINE CGRDH: Increase the leading dimension of H to ', &
              I0 )

!  end of subroutine CGRDH

      END SUBROUTINE CGRDH
