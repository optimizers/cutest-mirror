! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE UOFG ( data, n, X, f, G, GRAD )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n
      REAL ( KIND = wp ) :: f
      REAL ( KIND = wp ) :: X( n ), G( n )
      LOGICAL :: GRAD

!  Compute the value of the objective function and its gradient
!  for a function initially written in Standard Input Format (SIF).

!  G     is an array which gives the value of the gradient of
!        the objective function evaluated at X.
!        G(i) gives the partial derivative of the objective
!        function with respect to variable X(i).

!  Based on the subroutines ufn.f and ugr.f by Nick Gould, which are
!  in turn based on the subroutine SBMIN by Conn, Gould and Toint.

!  Ingrid Bongartz 
!  February 1993.


! ---------------------------------------------------------------------




! ---------------------------------------------------------------------



! ---------------------------------------------------------------------


! ---------------------------------------------------------------------

!  integer variables from the GLOBAL common block.

      INTEGER :: ESCALE, GSCALE, VSCALE, GVALS, XT, DGRAD 

!  integer variables from the LOCAL common block.


!  integer variables from the PRFCTS common block.


!  local variables.

      INTEGER :: i, j, ig, ifstat, igstat
      REAL ( KIND = wp ) :: ftt, DDOT, one, zero
      PARAMETER ( zero = 0.0_wp, one = 1.0_wp ) 
      EXTERNAL :: RANGE 

!  There are non-trivial group functions.

      DO 10 i = 1, MAX( data%nel, data%ng )
        data%ICALCF( i ) = i
   10 CONTINUE

!  Evaluate the element function values.

      CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), data%nel, &
                   data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                   data%IELVAR( 1 ), data%INTVAR( 1 ), &
                   data%ISTADH( 1 ), data%ISTEP( 1 ), &
                   data%ICALCF( 1 ),  &
                   data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                   data%lntvar, data%lintre, lfuval, data%lvscal, data%lepvlu,  &
                   1, ifstat )

!  Compute the group argument values ft.

      DO 100 ig = 1, data%ng
         ftt = - data%B( ig )

!  Include the contribution from the linear element 
!  only if the variable belongs to the first n variables.

         DO 30 i = data%ISTADA( ig ), data%ISTADA( ig + 1 ) - 1
            j = data%ICNA( i ) 
            IF ( j <= n )  &
               ftt = ftt + data%A( i ) * X( j )
   30    CONTINUE

!  Include the contributions from the nonlinear elements.

         DO 60 i = data%ISTADG( ig ), data%ISTADG( ig + 1 ) - 1
            ftt = ftt + &
                   data%ESCALE( i ) * data%FUVALS( data%IELING( i ) )
   60    CONTINUE
         data%FT( ig ) = ftt

!  Record the derivatives of trivial groups.

         IF ( data%GXEQX( ig ) ) data%GVALS( data%ng + ig ) = one
  100 CONTINUE

!  Compute the group function values.

!  All group functions are trivial.

      IF ( data%altriv ) THEN
!D       f = DDOT( data%ng, data%GSCALE( 1 ), 1, data%FT( 1 ), 1 )
      CALL DCOPY( data%ng, data%FT( 1 ), 1, data%GVALS( 1 ), 1 )
      CALL SETVL( data%ng, data%GVALS( data%ng + 1 ), 1, one )
      ELSE

!  Evaluate the group function values.

         CALL GROUP ( data%GVALS( 1 ), data%ng, data%FT( 1 ), &
                      data%GPVALU( 1 ), data%ng, &
                      data%ITYPEG( 1 ), data%ISTGP( 1 ), &
                      data%ICALCF( 1 ), &
                      data%lcalcg, data%ng1, data%lcalcg, data%lcalcg, data%lgpvlu, &
                      .FALSE., igstat )

!  Compute the objective function value.

         f = zero
         DO 110 ig = 1, data%ng
            IF ( data%GXEQX( ig ) ) THEN
               f = f + data%GSCALE( ig ) * data%FT( ig )
            ELSE
               f = f + data%GSCALE( ig ) * data%GVALS( ig )
            END IF
  110    CONTINUE
      END IF
      IF ( GRAD ) THEN

!  Evaluate the element function derivatives.

         CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), data%nel, &
                      data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                      data%IELVAR( 1 ), data%INTVAR( 1 ), &
                      data%ISTADH( 1 ), data%ISTEP( 1 ), &
                      data%ICALCF( 1 ),  &
                      data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                      data%lntvar, data%lintre, lfuval, data%lvscal, data%lepvlu,  &
                      2, ifstat )

!  Evaluate the group derivative values.

         IF ( .NOT. data%altriv ) CALL GROUP ( data%GVALS( 1 ), data%ng,  &
               data%FT( 1 ), data%GPVALU( 1 ), data%ng, &
               data%ITYPEG( 1 ), data%ISTGP( 1 ), &
               data%ICALCF( 1 ), &
               data%lcalcg, data%ng1, data%lcalcg, data%lcalcg, data%lgpvlu, &
               .TRUE., igstat )

!  Compute the gradient values. 

      CALL ELGRD( n, data%ng, data%firstg, data%ICNA( 1 ), data%licna, &
                      data%ISTADA( 1 ), data%lstada, data%IELING( 1 ), &
                      data%leling, data%ISTADG( 1 ), data%lstadg, &
                      data%ITYPEE( 1 ), data%lintre, &
                      data%ISTAEV( 1 ), data%lstaev, data%IELVAR( 1 ), &
                      data%lelvar, data%INTVAR( 1 ), data%lntvar, data%IWORK(LSVGRP + 1), &
                      data%lnvgrp, data%IWORK( data%lstajc + 1 ), data%lnstjc, &
                      data%IWORK( data%lstagv + 1 ), data%lnstgv, data%A( 1 ), data%la, &
                      data%GVALS( data%ng + 1 ), data%lgvals, &
                      data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ), &
                      data%GSCALE( 1 ), data%lgscal, &
                      data%ESCALE( 1 ), data%lescal, data%FUVALS( data%lgrjac + 1 ), &
                      data%lngrjc, data%WRK( 1 ), data%WRK( n + 1 ), data%maxsel, &
                      data%GXEQX( 1 ), data%lgxeqx, &
                      data%INTREP( 1 ), data%lintre, RANGE )
         data%firstg = .FALSE.

!  Store the gradient value.

         DO 300 i = 1, n
            G( i ) = data%FUVALS( data%lggfx + i )
  300    CONTINUE
      END IF

!  Update the counters for the report tool.

      data%nc2of = data%nc2of + 1
      IF( GRAD ) data%nc2og = data%nc2og + 1
      RETURN

!  end of UOFG.

      END
