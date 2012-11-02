! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE UFN ( data, N, X, F )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: N
      REAL ( KIND = wp ) :: F
      REAL ( KIND = wp ) :: X( N )

!  Compute the value of a groups partially separable function
!  initially written in Standard Input Format (SIF).

!  Based on the minimization subroutine SBMIN by Conn, Gould and Toint.

!  Nick Gould, for CGT productions.
!  July 1991.


! ---------------------------------------------------------------------




! ---------------------------------------------------------------------



! ---------------------------------------------------------------------


! ---------------------------------------------------------------------

!  integer variables from the GLOBAL common block.


!  integer variables from the LOCAL common block.




!  the common blocks


!  local variables.

      INTEGER :: I, J, IG, IFSTAT, IGSTAT
      REAL ( KIND = wp ) :: FTT, DDOT, ONE, ZERO
      PARAMETER ( ZERO = 0.0_wp, ONE = 1.0_wp )

!  increment the counter for calls to the objective function value

      data%nc2of = data%nc2of + 1

!  there are non-trivial group functions.

      DO 10 I = 1, MAX( data%nelnum, data%ng )
        data%ICALCF( I ) = I
   10 CONTINUE

!  evaluate the element function values.

      CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), data%nelnum, &
                   data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                   data%IELVAR( 1 ), data%INTVAR( 1 ), &
                   data%ISTADH( 1 ), data%ISTEP( 1 ), &
                   data%ICALCF( 1 ),  &
                   data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                   data%lntvar, data%lintre, LFUVAL, data%lvscal, data%lepvlu,  &
                   1, IFSTAT )

!  compute the group argument values ft.

      DO 100 IG = 1, data%ng
         FTT = - data%B( IG )

!  include the contribution from the linear element.

         DO 30 J = data%ISTADA( IG ), data%ISTADA( IG + 1 ) - 1
            FTT = FTT + data%A( J ) * X( data%ICNA( J ) )
   30    CONTINUE

!  include the contributions from the nonlinear elements.

         DO 60 J = data%ISTADG( IG ), data%ISTADG( IG + 1 ) - 1
            FTT = FTT + data%ESCALE( J ) * data%FUVALS( data%IELING( J ) )
   60    CONTINUE
         data%FT( IG ) = FTT
  100 CONTINUE

!  compute the group function values.

!  all group functions are trivial.

      IF ( data%altriv ) THEN
!D       F = DDOT( data%ng, data%GSCALE( 1 ), 1, data%FT( 1 ), 1 )
      CALL DCOPY( data%ng, data%FT( 1 ), 1, data%GVALS( 1 ), 1 )
      CALL DSETVL( data%ng, data%GVALS( data%ng + 1 ), 1, ONE )
      ELSE

!  evaluate the group function values.

         CALL GROUP ( data%GVALS( 1 ), data%ng, data%FT( 1 ), &
                      data%GPVALU( 1 ), data%ng, &
                      data%ITYPEG( 1 ), data%ISTGP( 1 ), &
                      data%ICALCF( 1 ), &
                      data%lcalcg, data%ng1, data%lcalcg, data%lcalcg, data%lgpvlu, &
                      .FALSE., IGSTAT )
         F = ZERO
         DO 220 IG = 1, data%ng
            IF ( data%GXEQX( IG ) ) THEN
               F = F + data%GSCALE( IG ) * data%FT( IG )
            ELSE
               F = F + data%GSCALE( IG ) * data%GVALS( IG )
            END IF
  220    CONTINUE
      END IF
      RETURN

!  end of UFN.

      END
