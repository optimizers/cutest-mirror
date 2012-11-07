! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE UFN ( data, n, X, f )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n
      REAL ( KIND = wp ) :: f
      REAL ( KIND = wp ) :: X( n )

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

      INTEGER :: i, j, ig, ifstat, igstat
      REAL ( KIND = wp ) :: ftt, DDOT, one, zero
      PARAMETER ( zero = 0.0_wp, one = 1.0_wp )

!  increment the counter for calls to the objective function value

      data%nc2of = data%nc2of + 1

!  there are non-trivial group functions.

      DO 10 i = 1, MAX( data%nel, data%ng )
        data%ICALCF( i ) = i
   10 CONTINUE

!  evaluate the element function values.

      CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), data%nel, &
                   data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                   data%IELVAR( 1 ), data%INTVAR( 1 ), &
                   data%ISTADH( 1 ), data%ISTEP( 1 ), &
                   data%ICALCF( 1 ),  &
                   data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                   data%lntvar, data%lintre, lfuval, data%lvscal, data%lepvlu,  &
                   1, ifstat )

!  compute the group argument values ft.

      DO 100 ig = 1, data%ng
         ftt = - data%B( ig )

!  include the contribution from the linear element.

         DO 30 j = data%ISTADA( ig ), data%ISTADA( ig + 1 ) - 1
            ftt = ftt + data%A( j ) * X( data%ICNA( j ) )
   30    CONTINUE

!  include the contributions from the nonlinear elements.

         DO 60 j = data%ISTADG( ig ), data%ISTADG( ig + 1 ) - 1
            ftt = ftt + data%ESCALE( j ) * data%FUVALS( data%IELING( j ) )
   60    CONTINUE
         data%FT( ig ) = ftt
  100 CONTINUE

!  compute the group function values.

!  all group functions are trivial.

      IF ( data%altriv ) THEN
!D       f = DDOT( data%ng, data%GSCALE( 1 ), 1, data%FT( 1 ), 1 )
      CALL DCOPY( data%ng, data%FT( 1 ), 1, data%GVALS( 1 ), 1 )
      CALL SETVL( data%ng, data%GVALS( data%ng + 1 ), 1, one )
      ELSE

!  evaluate the group function values.

         CALL GROUP ( data%GVALS( 1 ), data%ng, data%FT( 1 ), &
                      data%GPVALU( 1 ), data%ng, &
                      data%ITYPEG( 1 ), data%ISTGP( 1 ), &
                      data%ICALCF( 1 ), &
                      data%lcalcg, data%ng1, data%lcalcg, data%lcalcg, data%lgpvlu, &
                      .FALSE., igstat )
         f = zero
         DO 220 ig = 1, data%ng
            IF ( data%GXEQX( ig ) ) THEN
               f = f + data%GSCALE( ig ) * data%FT( ig )
            ELSE
               f = f + data%GSCALE( ig ) * data%GVALS( ig )
            END IF
  220    CONTINUE
      END IF
      RETURN

!  end of UFN.

      END
