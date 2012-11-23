! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CPROD( data, status, n, m, goth, X, Y, P, RESULT )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n, m
      INTEGER, INTENT( OUT ) :: status
      LOGICAL :: goth
      REAL ( KIND = wp ) :: X( n ), Y( m ), P( n ), RESULT( n )

!  Compute the matrix-vector product between the Hessian matrix
!  of the Lagrangian function for the problem and  a given vector P.
!  The result is placed in RESULT. If goth is .TRUE. the second
!  derivatives are assumed to have already been computed. If
!  the user is unsure, set goth = .FALSE. the first time a product
!  is required with the Hessian evaluated at X and V. X and V are
!  not used if goth = .TRUE.

!  Based on the minimization subroutine LANCELOT/SBMIN
!  by Conn, Gould and Toint.

!  Nick Gould, for CGT productions.
!  November, 1991.

!  Integer variables from the PRFCTS common block.


!  Local variables

      INTEGER :: i, ig, j, ifstat, igstat
      REAL ( KIND = wp ) :: ftt
      EXTERNAL :: RANGE 

!  there are non-trivial group functions

      IF ( .NOT. goth ) THEN
        DO i = 1, MAX( data%nel, data%ng )
          data%ICALCF( i ) = i
        END DO

!  evaluate the element function values

        CALL ELFUN( data%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,        &
                    data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,        &
                    data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,         &
                    data%lelvar, data%lntvar, data%lstadh, data%lstep,         &
                    data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,        &
                    1, ifstat )
        IF ( ifstat /= 0 ) GO TO 930

!  evaluate the element function gradient and Hessian values

        CALL ELFUN( data%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,        &
                    data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,        &
                    data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,         &
                    data%lelvar, data%lntvar, data%lstadh, data%lstep,         &
                    data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,        &
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
          CALL GROUP( data%GVALS, data%ng, data%FT, data%GPVALU, data%ng,      &
                      data%ITYPEG, data%ISTGP, data%ICALCF, data%ltypeg,       &
                      data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,       &
                      .TRUE., igstat )
          IF ( igstat /= 0 ) GO TO 930
        END IF

!  define the real work space needed for ELGRD. Ensure that there is 
!  sufficient space

        IF ( data%lwk2 < data%ng ) THEN
          IF ( data%out > 0 ) WRITE( data%out, 2000 )
          status = 2 ; RETURN
        END IF

!  change the group weightings to include the contributions from the
!  Lagrange multipliers.

        IF ( data%numcon > 0 ) THEN
          DO ig = 1, data%ng
            i = data%KNDOFC( ig )
            IF ( i == 0 ) THEN
              data%GSCALE_used( ig ) = data%GSCALE( ig )
            ELSE
              data%GSCALE_used( ig ) = data%GSCALE( ig ) * Y( i )
            END IF
          END DO

!  compute the gradient value

          CALL CUTEST_form_gradients( n, data%ng, data%nel, data%ntotel,       &
                 data%nvrels, data%nnza, data%nvargp, data%firstg, data%ICNA,  &
                 data%ISTADA, data%IELING, data%ISTADG, data%ISTAEV,           &
                 data%IELVAR, data%INTVAR, data%A, data%GVALS( : , 2 ),        &
                 data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ),      &
                 data%GSCALE_used, data%ESCALE, data%FUVALS( data%lgrjac + 1 ),&
                 data%GXEQX, data%INTREP, data%ISVGRP, data%ISTAGV,            &
                 data%ITYPEE, data%ISTAJC, data%W_ws, data%W_el, RANGE )
        ELSE
          CALL CUTEST_form_gradients( n, data%ng, data%nel, data%ntotel,       &
                 data%nvrels, data%nnza, data%nvargp, data%firstg, data%ICNA,  &
                 data%ISTADA, data%IELING, data%ISTADG, data%ISTAEV,           &
                 data%IELVAR, data%INTVAR, data%A, data%GVALS( : , 2 ),        &
                 data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ),      &
                 data%GSCALE, data%ESCALE, data%FUVALS( data%lgrjac + 1 ),     &
                 data%GXEQX, data%INTREP, data%ISVGRP, data%ISTAGV,            &
                 data%ITYPEE, data%ISTAJC, data%W_ws, data%W_el, RANGE )
        END IF

!  Compute the gradient value

!          CALL ELGRD( n, data%ng, data%firstg, data%ICNA( 1 ), data%licna, &
!                         data%ISTADA( 1 ), data%lstada, data%IELING( 1 ), &
!                         data%leling, data%ISTADG( 1 ), data%lstadg, &
!                         data%ITYPEE( 1 ), data%lintre, &
!                         data%ISTAEV( 1 ), data%lstaev, data%IELVAR( 1 ), &
!                         data%lelvar, data%INTVAR( 1 ), data%lntvar, &
!                         data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, data%IWORK( data%lstajc + 1 ), &
!                         data%lnstjc, data%IWORK( data%lstagv + 1 ), data%lnstgv, &
!                         data%A( 1 ), data%la, data%GVALS( : , 2 ), data%lgvals, &
!                         data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ), &
!                         data%WRK( 1 ), data%ng, &
!                         data%ESCALE( 1 ), data%lescal, data%FUVALS( data%lgrjac + 1 ), &
!                         data%lngrjc, data%WRK( 1 ), data%WRK( n + 1 ), &
!                         data%maxsel, data%GXEQX( 1 ), data%lgxeqx, &
!                         data%INTREP( 1 ), data%lintre, RANGE )
!        ELSE
!          CALL ELGRD( n, data%ng, data%firstg, data%ICNA( 1 ), data%licna, &
!                         data%ISTADA( 1 ), data%lstada, data%IELING( 1 ), &
!                         data%leling, data%ISTADG( 1 ), data%lstadg, &
!                         data%ITYPEE( 1 ), data%lintre, &
!                         data%ISTAEV( 1 ), data%lstaev, data%IELVAR( 1 ), &
!                         data%lelvar, data%INTVAR( 1 ), data%lntvar, &
!                         data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, data%IWORK( data%lstajc + 1 ), &
!                         data%lnstjc, data%IWORK( data%lstagv + 1 ), data%lnstgv, &
!                         data%A( 1 ), data%la, data%GVALS( : , 2 ), data%lgvals, &
!                         data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ), &
!                         data%GSCALE( 1 ), data%lgscal, &
!                         data%ESCALE( 1 ), data%lescal, data%FUVALS( data%lgrjac + 1 ), &
!                         data%lngrjc, data%WRK( 1 ), data%WRK( n + 1 ), &
!                         data%maxsel, data%GXEQX( 1 ), data%lgxeqx, &
!                         data%INTREP( 1 ), data%lintre, RANGE )
        data%firstg = .FALSE.
      END IF

!  ensure that the product involves all components of P

      DO i = 1, n
        data%IVAR( i ) = i
!       data%IWORK( data%lnnonz + i ) = i
      END DO

!  initialize RESULT as the zero vector

      RESULT( : n ) = 0.0_wp

!  define the real work space needed for HSPRD.  Ensure that there is 
!  sufficient space

!     nn = data%ninvar + n
!     lnwk = MAX( data%ng, data%maxsel )
!     lnwkb = data%maxsin
!     lnwkc = data%maxsin
!     lwkb = lnwk
!     lwkc = lwkb + lnwkb

!  Evaluate the product.

      IF ( data%numcon > 0 ) THEN
        CALL CUTEST_hessian_times_vector(                                      &
          data%n, data%ng, data%nel, data%ntotel, data%nvrels, data%nvargp,    &
          n, 1, n, data%nnonnz, data%nbprod, data%alllin,                      &
          data%IVAR, data%ISTAEV, data%ISTADH, data%INTVAR, data%IELING,       &
          data%IELVAR, data%ISWKSP( : data%ntotel ), data%INNONZ( : n ),       &
          data%P, data%Q, data%GVALS( : , 2 ) , data%GVALS( : , 3 ),           &
          data%GRJAC, data%GSCALE_used,                                        &
          data%ESCALE, data%FUVALS( : data%lnhuvl ),                           &
          data%lnhuvl, data%GXEQX, data%INTREP, .TRUE., data%IGCOLJ,           &
          data%ISLGRP, data%ISVGRP, data%ISTAGV, data%IVALJR, data%ITYPEE,     &
          data%ISYMMH, data%ISTAJC, data%IUSED, data%LIST_elements,            &
          data%LINK_elem_uses_var, data%NZ_comp_w, data%W_ws, data%W_el,       &
          data%W_in, data%H_in, RANGE, data%skipg, data%KNDOFC )
      ELSE
        CALL CUTEST_hessian_times_vector(                                      &
          data%n, data%ng, data%nel, data%ntotel, data%nvrels, data%nvargp,    &
          n, 1, n, data%nnonnz, data%nbprod, data%alllin,                      &
          data%IVAR, data%ISTAEV, data%ISTADH, data%INTVAR, data%IELING,       &
          data%IELVAR, data%ISWKSP( : data%ntotel ), data%INNONZ( : n ),       &
          data%P, data%Q, data%GVALS( : , 2 ) , data%GVALS( : , 3 ),           &
          data%GRJAC, data%GSCALE,                                             &
          data%ESCALE, data%FUVALS( : data%lnhuvl ),                           &
          data%lnhuvl, data%GXEQX, data%INTREP, .TRUE., data%IGCOLJ,           &
          data%ISLGRP, data%ISVGRP, data%ISTAGV, data%IVALJR, data%ITYPEE,     &
          data%ISYMMH, data%ISTAJC, data%IUSED, data%LIST_elements,            &
          data%LINK_elem_uses_var, data%NZ_comp_w, data%W_ws, data%W_el,       &
          data%W_in, data%H_in, RANGE, data%skipg, data%KNDOFC )
      END IF

!      IF ( data%numcon > 0 ) THEN
!        CALL DHSPRD( n, nn, data%ng, data%ntotel, n, 1, n, nbprod, 
!               data%nel == 0, &
!               data%IVAR( 1 ), data%ISTAEV( 1 ), data%lstaev, &
!               data%ISTADH( 1 ), data%lstadh, data%INTVAR( 1 ), &
!               data%lntvar, data%IELING( 1 ), data%leling, 
!               data%IELVAR( 1 ), &
!               data%lelvar, data%IWORK( data%lstajc + 1 ), 
!               data%lnstjc, data%IWORK( data%lselts + 1 ), &
!               data%lnelts, data%IWORK( data%lsptrs + 1 ), 
!               data%lnptrs, data%IWORK( data%lgcolj + 1 ),  &
!               data%lngclj, data%IWORK( data%lslgrp + 1 ), 
!               data%lnlgrp, data%IWORK( data%lswksp + 1 ),  &
!               data%lnwksp, data%IWORK( data%lsvgrp + 1 ), 
!               data%lnvgrp, data%IWORK( data%lstagv + 1 ), &
!               data%lnstgv, data%IWORK( data%lvaljr + 1 ), 
!               data%lnvljr, data%ITYPEE( 1 ),  &
!               data%lintre, nnonnz, data%IWORK( data%lnnonz + 1 ), 
!               data%lnnnon, &
!               data%IWORK( data%liused + 1 ), data%lniuse, 
!               data%IWORK( data%lnonz2 + 1 ), &
!               data%lnnno2, data%IWORK( data%lsymmh + 1 ), 
!               data%maxsin, P, RESULT, &
!               data%GVALS( : , 2 ), data%GVALS( : , 3 ), &
!               data%FUVALS( data%lgrjac + 1 ), data%lngrjc, data%WRK( 1 ), &
!               data%ESCALE( 1 ), data%lescal, data%FUVALS, data%lnhuvl, &
!               data%WRK( 1 ), lnwk, data%WRK( lwkb + 1 ), &
!               lnwkb, data%WRK( lwkc + 1 ), lnwkc, &
!               data%GXEQX( 1 ), data%lgxeqx, data%INTREP( 1 ), &
!               data%lintre, .TRUE., RANGE )
!      ELSE
!        CALL DHSPRD( n, nn, data%ng, data%ntotel, n, 1, n, nbprod, 
!               data%nel == 0, &
!               data%IVAR( 1 ), data%ISTAEV( 1 ), data%lstaev, &
!               data%ISTADH( 1 ), data%lstadh, data%INTVAR( 1 ), &
!               data%lntvar, data%IELING( 1 ), data%leling, 
!               data%IELVAR( 1 ), &
!               data%lelvar, data%IWORK( data%lstajc + 1 ), 
!               data%lnstjc, data%IWORK( data%lselts + 1 ), &
!               data%lnelts, data%IWORK( data%lsptrs + 1 ), 
!               data%lnptrs, data%IWORK( data%lgcolj + 1 ),  &
!               data%lngclj, data%IWORK( data%lslgrp + 1 ), 
!               data%lnlgrp, data%IWORK( data%lswksp + 1 ),  &
!               data%lnwksp, data%IWORK( data%lsvgrp + 1 ), 
!               data%lnvgrp, data%IWORK( data%lstagv + 1 ), &
!               data%lnstgv, data%IWORK( data%lvaljr + 1 ), 
!               data%lnvljr, data%ITYPEE( 1 ),  &
!               data%lintre, nnonnz, data%IWORK( data%lnnonz + 1 ), 
!               data%lnnnon, &
!               data%IWORK( data%liused + 1 ), data%lniuse, 
!               data%IWORK( data%lnonz2 + 1 ), &
!               data%lnnno2, data%IWORK( data%lsymmh + 1 ), 
!               data%maxsin, P, RESULT, &
!               data%GVALS( : , 2 ), data%GVALS( : , 3 ), &
!               data%FUVALS( data%lgrjac + 1 ), data%lngrjc, 
!               data%GSCALE( 1 ), &
!               data%ESCALE( 1 ), data%lescal, data%FUVALS, data%lnhuvl, &
!               data%WRK( 1 ), lnwk, data%WRK( lwkb + 1 ), &
!               lnwkb, data%WRK( lwkc + 1 ), lnwkc, &
!               data%GXEQX( 1 ), data%lgxeqx, data%INTREP( 1 ), &
!               data%lintre, .TRUE., RANGE )
!      END IF

!  update the counters for the report tool

      data%nhvpr = data%nhvpr + 1
      IF ( .NOT. goth ) THEN
        data%nc2oh = data%nc2oh + 1
        data%nc2ch = data%nc2ch + data%pnc
      END IF
      status = 0
      RETURN

!  unsuccessful returns

  930 CONTINUE
      IF ( data%out > 0 ) WRITE( data%out,                                     &
        "( ' ** SUBROUTINE CPROD: error flag raised during SIF evaluation' )" )
      status = 3
      RETURN

!  non-executable statements

 2000 FORMAT( ' ** SUBROUTINE CPROD: Increase the size of WK' )

!  end of subroutine CPROD

      END SUBROUTINE CPROD


