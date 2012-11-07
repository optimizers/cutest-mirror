! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CPROD( data, n, m, GOTH, X, lv, V, P, RESULT )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n, m, lv
      LOGICAL :: GOTH
      REAL ( KIND = wp ) :: X( n ), V ( lv ), P( n ), RESULT( n )

!  Compute the matrix-vector product between the Hessian matrix
!  of the Lagrangian function for the problem and  a given vector P.
!  The result is placed in RESULT. If GOTH is .TRUE. the second
!  derivatives are assumed to have already been computed. If
!  the user is unsure, set GOTH = .FALSE. the first time a product
!  is required with the Hessian evaluated at X and V. X and V are
!  not used if GOTH = .TRUE.

!  Based on the minimization subroutine data%laNCELOT/SBMIN
!  by Conn, Gould and Toint.

!  Nick Gould, for CGT productions.
!  November, 1991.

!  Integer variables from the PRFCTS common block.


!  Local variables

      INTEGER :: i, ig, j, nn, nbprod, nnonnz
      INTEGER :: lnwk, lnwkb, lnwkc, lwkb, lwkc
      INTEGER :: ifstat, igstat
      REAL ( KIND = wp ) :: zero, one, ftt
      PARAMETER ( zero = 0.0_wp, one = 1.0_wp )
      EXTERNAL :: RANGE 

!  There are non-trivial group functions.

      IF ( .NOT. GOTH ) THEN
         DO 10 i = 1, MAX( data%nel, data%ng )
           data%ICALCF( i ) = i
   10    CONTINUE

!  Evaluate the element function values.

         CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), data%nel, &
                      data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                      data%IELVAR( 1 ), data%INTVAR( 1 ), &
                      data%ISTADH( 1 ), data%ISTEP( 1 ), &
                      data%ICALCF( 1 ),  &
                      data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                      data%lntvar, data%lintre, lfuval, data%lvscal, data%lepvlu,  &
                      1, ifstat )

!  Evaluate the element function values.

         CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), data%nel, &
                      data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                      data%IELVAR( 1 ), data%INTVAR( 1 ), &
                      data%ISTADH( 1 ), data%ISTEP( 1 ), &
                      data%ICALCF( 1 ),  &
                      data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                      data%lntvar, data%lintre, lfuval, data%lvscal, data%lepvlu,  &
                      3, ifstat )

!  Compute the group argument values ft.

         DO 70 ig = 1, data%ng
            ftt = - data%B( ig )

!  Include the contribution from the linear element.

            DO 30 j = data%ISTADA( ig ), data%ISTADA( ig + 1 ) - 1
               ftt = ftt + data%A( j ) * X( data%ICNA( j ) )
   30       CONTINUE

!  Include the contributions from the nonlinear elements.

            DO 60 j = data%ISTADG( ig ), data%ISTADG( ig + 1 ) - 1
               ftt = ftt + data%ESCALE( j ) * data%FUVALS( data%IELING( J))
   60       CONTINUE
            data%FT( ig ) = ftt

!  Record the derivatives of trivial groups.

            IF ( data%GXEQX( ig ) ) THEN
               data%GVALS( data%ng + ig ) = one
               data%GVALS( 2 * data%ng + ig ) = zero
            END IF
   70    CONTINUE

!  Evaluate the group derivative values.

         IF ( .NOT. data%altriv ) CALL GROUP ( data%GVALS( 1 ), data%ng, &
               data%FT( 1 ), data%GPVALU( 1 ), data%ng, &
               data%ITYPEG( 1 ), data%ISTGP( 1 ), &
               data%ICALCF( 1 ), &
               data%lcalcg, data%ng1, data%lcalcg, data%lcalcg, data%lgpvlu, &
               .TRUE., igstat )

!  Define the real work space needed for ELGRD.
!  Ensure that there is sufficient space.

         IF ( data%lwk2 < data%ng ) THEN
            IF ( iout > 0 ) WRITE( iout, 2000 )
            STOP
         END IF
!                            for unconstrained problems.
         IF ( data%numcon > 0 ) THEN

!  Change the group weightings to include the contributions from
!  the Lagrange multipliers.

            DO 80 ig = 1, data%ng
               i = data%KNDOFC( ig )
               IF ( i == 0 ) THEN
                  data%WRK( ig ) = data%GSCALE( ig )
               ELSE
                  data%WRK( ig ) = data%GSCALE( ig ) * V( i )
               END IF
   80       CONTINUE

!  Compute the gradient value.

      CALL ELGRD( n, data%ng, data%firstg, data%ICNA( 1 ), data%licna, &
                         data%ISTADA( 1 ), data%lstada, data%IELING( 1 ), &
                         data%leling, data%ISTADG( 1 ), data%lstadg, &
                         data%ITYPEE( 1 ), data%lintre, &
                         data%ISTAEV( 1 ), data%lstaev, data%IELVAR( 1 ), &
                         data%lelvar, data%INTVAR( 1 ), data%lntvar, &
                         data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, data%IWORK( data%lstajc + 1 ), &
                         data%lnstjc, data%IWORK( data%lstagv + 1 ), data%lnstgv, &
                         data%A( 1 ), data%la, data%GVALS( data%ng + 1 ), data%lgvals, &
                         data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ), &
                         data%WRK( 1 ), data%ng, &
                         data%ESCALE( 1 ), data%lescal, data%FUVALS( data%lgrjac + 1 ), &
                         data%lngrjc, data%WRK( 1 ), data%WRK( n + 1 ), &
                         data%maxsel, data%GXEQX( 1 ), data%lgxeqx, &
                         data%INTREP( 1 ), data%lintre, RANGE )
         ELSE

!  Compute the gradient value.

      CALL ELGRD( n, data%ng, data%firstg, data%ICNA( 1 ), data%licna, &
                         data%ISTADA( 1 ), data%lstada, data%IELING( 1 ), &
                         data%leling, data%ISTADG( 1 ), data%lstadg, &
                         data%ITYPEE( 1 ), data%lintre, &
                         data%ISTAEV( 1 ), data%lstaev, data%IELVAR( 1 ), &
                         data%lelvar, data%INTVAR( 1 ), data%lntvar, &
                         data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, data%IWORK( data%lstajc + 1 ), &
                         data%lnstjc, data%IWORK( data%lstagv + 1 ), data%lnstgv, &
                         data%A( 1 ), data%la, data%GVALS( data%ng + 1 ), data%lgvals, &
                         data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ), &
                         data%GSCALE( 1 ), data%lgscal, &
                         data%ESCALE( 1 ), data%lescal, data%FUVALS( data%lgrjac + 1 ), &
                         data%lngrjc, data%WRK( 1 ), data%WRK( n + 1 ), &
                         data%maxsel, data%GXEQX( 1 ), data%lgxeqx, &
                         data%INTREP( 1 ), data%lintre, RANGE )
         END IF
         data%firstg = .FALSE.
      END IF

!  Ensure that the product involves all components of P.

      DO 100 i = 1, n
         data%IVAR( i ) = i
         data%IWORK( data%lnnonz + i ) = i
  100 CONTINUE

!  Initialize RESULT as the zero vector.

      CALL SETVL( n, RESULT, 1, zero )

!  Define the real work space needed for HSPRD.
!  Ensure that there is sufficient space.

      nn = data%ninvar + n
      lnwk = MAX( data%ng, data%maxsel )
      lnwkb = data%maxsin
      lnwkc = data%maxsin
      lwkb = lnwk
      lwkc = lwkb + lnwkb

!  Evaluate the product.


!            CALL HSPRD_hessian_times_vector(                                  &
!                n , ng, nel, S%ntotel, S%nvrels, S%nvargp,                    &
!                inform%nvar  , nvar1 , S%nvar2 , S%nnonnz,                    &
!                S%nbprod, S%alllin, IVAR , ISTAEV, ISTADH, INTVAR, IELING,    &
!                IELVAR, ISWKSP( : S%ntotel ), INNONZ( : n ),                  &
!                P , Q , GVALS( : , 2 )  , GVALS( : , 3 ),                     &
!                GRJAC, GSCALE_used, ESCALE, FUVALS( : S%lnhuvl ), S%lnhuvl,   &
!                GXEQX_used , INTREP, S%densep,                                &
!                IGCOLJ, ISLGRP, ISVGRP, ISTAGV, IVALJR, ITYPEE, ISYMMH,       &
!                ISTAJC, IUSED, LIST_elements, LINK_elem_uses_var,             &
!                NZ_comp_w, W_ws, W_el, W_in, H_in, RANGE, S%skipg, KNDOFG )




      IF ( data%numcon > 0 ) THEN
      CALL DHSPRD( n, nn, data%ng, data%ntotel, n, 1, n, nbprod, data%nel == 0, &
                   data%IVAR( 1 ), data%ISTAEV( 1 ), data%lstaev, &
                   data%ISTADH( 1 ), data%lstadh, data%INTVAR( 1 ), &
                   data%lntvar, data%IELING( 1 ), data%leling, data%IELVAR( 1 ), &
                   data%lelvar, data%IWORK( data%lstajc + 1 ), data%lnstjc, data%IWORK( data%lselts + 1 ), &
                   data%lnelts, data%IWORK( data%lsptrs + 1 ), data%lnptrs, data%IWORK( data%lgcolj + 1 ),  &
                   data%lngclj, data%IWORK( data%lslgrp + 1 ), data%lnlgrp, data%IWORK( data%lswksp + 1 ),  &
                   data%lnwksp, data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, data%IWORK( data%lstagv + 1 ), &
                   data%lnstgv, data%IWORK( data%lvaljr + 1 ), data%lnvljr, data%ITYPEE( 1 ),  &
                   data%lintre, nnonnz, data%IWORK( data%lnnonz + 1 ), data%lnnnon, &
                   data%IWORK( data%liused + 1 ), data%lniuse, data%IWORK( data%lnonz2 + 1 ), &
                   data%lnnno2, data%IWORK( data%lsymmh + 1 ), data%maxsin, P, RESULT, &
                   data%GVALS( data%ng + 1 ), data%GVALS( 2 * data%ng + 1 ), &
                   data%FUVALS( data%lgrjac + 1 ), data%lngrjc, data%WRK( 1 ), &
                   data%ESCALE( 1 ), data%lescal, data%FUVALS, data%lnhuvl, &
                   data%WRK( 1 ), lnwk, data%WRK( lwkb + 1 ), &
                   lnwkb, data%WRK( lwkc + 1 ), lnwkc, &
                   data%GXEQX( 1 ), data%lgxeqx, data%INTREP( 1 ), &
                   data%lintre, .TRUE., RANGE )
      ELSE
      CALL DHSPRD( n, nn, data%ng, data%ntotel, n, 1, n, nbprod, data%nel == 0, &
                   data%IVAR( 1 ), data%ISTAEV( 1 ), data%lstaev, &
                   data%ISTADH( 1 ), data%lstadh, data%INTVAR( 1 ), &
                   data%lntvar, data%IELING( 1 ), data%leling, data%IELVAR( 1 ), &
                   data%lelvar, data%IWORK( data%lstajc + 1 ), data%lnstjc, data%IWORK( data%lselts + 1 ), &
                   data%lnelts, data%IWORK( data%lsptrs + 1 ), data%lnptrs, data%IWORK( data%lgcolj + 1 ),  &
                   data%lngclj, data%IWORK( data%lslgrp + 1 ), data%lnlgrp, data%IWORK( data%lswksp + 1 ),  &
                   data%lnwksp, data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, data%IWORK( data%lstagv + 1 ), &
                   data%lnstgv, data%IWORK( data%lvaljr + 1 ), data%lnvljr, data%ITYPEE( 1 ),  &
                   data%lintre, nnonnz, data%IWORK( data%lnnonz + 1 ), data%lnnnon, &
                   data%IWORK( data%liused + 1 ), data%lniuse, data%IWORK( data%lnonz2 + 1 ), &
                   data%lnnno2, data%IWORK( data%lsymmh + 1 ), data%maxsin, P, RESULT, &
                   data%GVALS( data%ng + 1 ), data%GVALS( 2 * data%ng + 1 ), &
                   data%FUVALS( data%lgrjac + 1 ), data%lngrjc, data%GSCALE( 1 ), &
                   data%ESCALE( 1 ), data%lescal, data%FUVALS, data%lnhuvl, &
                   data%WRK( 1 ), lnwk, data%WRK( lwkb + 1 ), &
                   lnwkb, data%WRK( lwkc + 1 ), lnwkc, &
                   data%GXEQX( 1 ), data%lgxeqx, data%INTREP( 1 ), &
                   data%lintre, .TRUE., RANGE )
      END IF

!  Update the counters for the report tool.

      data%nhvpr = data%nhvpr + 1
      IF ( .NOT. GOTH ) THEN
         data%nc2oh = data%nc2oh + 1
         data%nc2ch = data%nc2ch + data%pnc
      END IF
      RETURN

! Non-executable statements.

 2000 FORMAT( ' ** SUBROUTINE CPROD: Increase the size of WK ' )

!  end of CPROD.

      END


