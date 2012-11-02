! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CPROD( data, N, M, GOTH, X, LV, V, P, RESULT )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: N, M, LV
      LOGICAL :: GOTH
      REAL ( KIND = wp ) :: X( N ), V ( LV ), P( N ), RESULT( N )

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


! ---------------------------------------------------------------------




! ---------------------------------------------------------------------



! ---------------------------------------------------------------------


! ---------------------------------------------------------------------





!  Integer variables from the PRFCTS common block.


!  Local variables

      INTEGER :: I, IG, J, NN, NBPROD, NNONNZ
      INTEGER :: LNWK, LNWKB, LNWKC, LWKB, LWKC
      INTEGER :: IFSTAT, IGSTAT
      REAL ( KIND = wp ) :: ZERO, ONE, FTT
      PARAMETER ( ZERO = 0.0_wp, ONE = 1.0_wp )
      EXTERNAL :: RANGE 

!  There are non-trivial group functions.

      IF ( .NOT. GOTH ) THEN
         DO 10 I = 1, MAX( data%nelnum, data%ng )
           data%ICALCF( I ) = I
   10    CONTINUE

!  Evaluate the element function values.

         CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), data%nelnum, &
                      data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                      data%IELVAR( 1 ), data%INTVAR( 1 ), &
                      data%ISTADH( 1 ), data%ISTEP( 1 ), &
                      data%ICALCF( 1 ),  &
                      data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                      data%lntvar, data%lintre, LFUVAL, data%lvscal, data%lepvlu,  &
                      1, IFSTAT )

!  Evaluate the element function values.

         CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), data%nelnum, &
                      data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                      data%IELVAR( 1 ), data%INTVAR( 1 ), &
                      data%ISTADH( 1 ), data%ISTEP( 1 ), &
                      data%ICALCF( 1 ),  &
                      data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                      data%lntvar, data%lintre, LFUVAL, data%lvscal, data%lepvlu,  &
                      3, IFSTAT )

!  Compute the group argument values ft.

         DO 70 IG = 1, data%ng
            FTT = - data%B( IG )

!  Include the contribution from the linear element.

            DO 30 J = data%ISTADA( IG ), data%ISTADA( IG + 1 ) - 1
               FTT = FTT + data%A( J ) * X( data%ICNA( J ) )
   30       CONTINUE

!  Include the contributions from the nonlinear elements.

            DO 60 J = data%ISTADG( IG ), data%ISTADG( IG + 1 ) - 1
               FTT = FTT + data%ESCALE( J ) * data%FUVALS( data%IELING( J))
   60       CONTINUE
            data%FT( IG ) = FTT

!  Record the derivatives of trivial groups.

            IF ( data%GXEQX( IG ) ) THEN
               data%GVALS( data%ng + IG ) = ONE
               data%GVALS( 2 * data%ng + IG ) = ZERO
            END IF
   70    CONTINUE

!  Evaluate the group derivative values.

         IF ( .NOT. data%altriv ) CALL GROUP ( data%GVALS( 1 ), data%ng, &
               data%FT( 1 ), data%GPVALU( 1 ), data%ng, &
               data%ITYPEG( 1 ), data%ISTGP( 1 ), &
               data%ICALCF( 1 ), &
               data%lcalcg, data%ng1, data%lcalcg, data%lcalcg, data%lgpvlu, &
               .TRUE., IGSTAT )

!  Define the real work space needed for ELGRD.
!  Ensure that there is sufficient space.

         IF ( data%lwk2 < data%ng ) THEN
            IF ( IOUT > 0 ) WRITE( IOUT, 2000 )
            STOP
         END IF
!                            for unconstrained problems.
         IF ( data%numcon > 0 ) THEN

!  Change the group weightings to include the contributions from
!  the Lagrange multipliers.

            DO 80 IG = 1, data%ng
               I = data%KNDOFC( IG )
               IF ( I == 0 ) THEN
                  data%WRK( IG ) = data%GSCALE( IG )
               ELSE
                  data%WRK( IG ) = data%GSCALE( IG ) * V( I )
               END IF
   80       CONTINUE

!  Compute the gradient value.

      CALL DELGRD( N, data%ng, data%firstg, data%ICNA( 1 ), data%licna, &
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
                         data%lngrjc, data%WRK( 1 ), data%WRK( N + 1 ), &
                         data%maxsel, data%GXEQX( 1 ), data%lgxeqx, &
                         data%INTREP( 1 ), data%lintre, RANGE )
         ELSE

!  Compute the gradient value.

      CALL DELGRD( N, data%ng, data%firstg, data%ICNA( 1 ), data%licna, &
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
                         data%lngrjc, data%WRK( 1 ), data%WRK( N + 1 ), &
                         data%maxsel, data%GXEQX( 1 ), data%lgxeqx, &
                         data%INTREP( 1 ), data%lintre, RANGE )
         END IF
         data%firstg = .FALSE.
      END IF

!  Ensure that the product involves all components of P.

      DO 100 I = 1, N
         data%IVAR( I ) = I
         data%IWORK( data%lnnonz + I ) = I
  100 CONTINUE

!  Initialize RESULT as the zero vector.

      CALL DSETVL( N, RESULT, 1, ZERO )

!  Define the real work space needed for HSPRD.
!  Ensure that there is sufficient space.

      NN = data%ninvar + N
      LNWK = MAX( data%ng, data%maxsel )
      LNWKB = data%maxsin
      LNWKC = data%maxsin
      LWKB = LNWK
      LWKC = LWKB + LNWKB

!  Evaluate the product.

      IF ( data%numcon > 0 ) THEN
      CALL DHSPRD( N, NN, data%ng, data%ngel, N, 1, N, NBPROD, data%nelnum == 0, &
                   data%IVAR( 1 ), data%ISTAEV( 1 ), data%lstaev, &
                   data%ISTADH( 1 ), data%lstadh, data%INTVAR( 1 ), &
                   data%lntvar, data%IELING( 1 ), data%leling, data%IELVAR( 1 ), &
                   data%lelvar, data%IWORK( data%lstajc + 1 ), data%lnstjc, data%IWORK( data%lselts + 1 ), &
                   data%lnelts, data%IWORK( data%lsptrs + 1 ), data%lnptrs, data%IWORK( data%lgcolj + 1 ),  &
                   data%lngclj, data%IWORK( data%lslgrp + 1 ), data%lnlgrp, data%IWORK( data%lswksp + 1 ),  &
                   data%lnwksp, data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, data%IWORK( data%lstagv + 1 ), &
                   data%lnstgv, data%IWORK( data%lvaljr + 1 ), data%lnvljr, data%ITYPEE( 1 ),  &
                   data%lintre, NNONNZ, data%IWORK( data%lnnonz + 1 ), data%lnnnon, &
                   data%IWORK( data%liused + 1 ), data%lniuse, data%IWORK( data%lnonz2 + 1 ), &
                   data%lnnno2, data%IWORK( data%lsymmh + 1 ), data%maxsin, P, RESULT, &
                   data%GVALS( data%ng + 1 ), data%GVALS( 2 * data%ng + 1 ), &
                   data%FUVALS( data%lgrjac + 1 ), data%lngrjc, data%WRK( 1 ), &
                   data%ESCALE( 1 ), data%lescal, data%FUVALS, data%lnhuvl, &
                   data%WRK( 1 ), LNWK, data%WRK( LWKB + 1 ), &
                   LNWKB, data%WRK( LWKC + 1 ), LNWKC, &
                   data%GXEQX( 1 ), data%lgxeqx, data%INTREP( 1 ), &
                   data%lintre, .TRUE., RANGE )
      ELSE
      CALL DHSPRD( N, NN, data%ng, data%ngel, N, 1, N, NBPROD, data%nelnum == 0, &
                   data%IVAR( 1 ), data%ISTAEV( 1 ), data%lstaev, &
                   data%ISTADH( 1 ), data%lstadh, data%INTVAR( 1 ), &
                   data%lntvar, data%IELING( 1 ), data%leling, data%IELVAR( 1 ), &
                   data%lelvar, data%IWORK( data%lstajc + 1 ), data%lnstjc, data%IWORK( data%lselts + 1 ), &
                   data%lnelts, data%IWORK( data%lsptrs + 1 ), data%lnptrs, data%IWORK( data%lgcolj + 1 ),  &
                   data%lngclj, data%IWORK( data%lslgrp + 1 ), data%lnlgrp, data%IWORK( data%lswksp + 1 ),  &
                   data%lnwksp, data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, data%IWORK( data%lstagv + 1 ), &
                   data%lnstgv, data%IWORK( data%lvaljr + 1 ), data%lnvljr, data%ITYPEE( 1 ),  &
                   data%lintre, NNONNZ, data%IWORK( data%lnnonz + 1 ), data%lnnnon, &
                   data%IWORK( data%liused + 1 ), data%lniuse, data%IWORK( data%lnonz2 + 1 ), &
                   data%lnnno2, data%IWORK( data%lsymmh + 1 ), data%maxsin, P, RESULT, &
                   data%GVALS( data%ng + 1 ), data%GVALS( 2 * data%ng + 1 ), &
                   data%FUVALS( data%lgrjac + 1 ), data%lngrjc, data%GSCALE( 1 ), &
                   data%ESCALE( 1 ), data%lescal, data%FUVALS, data%lnhuvl, &
                   data%WRK( 1 ), LNWK, data%WRK( LWKB + 1 ), &
                   LNWKB, data%WRK( LWKC + 1 ), LNWKC, &
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


