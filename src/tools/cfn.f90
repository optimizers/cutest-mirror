! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CFN ( data, N, M, X, F, LC, C )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: N, M, LC
      REAL ( KIND = wp ) :: F
      REAL ( KIND = wp ) :: X ( N ), C ( LC )

!  Compute the values of the objective function and general constraints
!  of a function initially written in Standard Input Format (SIF).

!  Based on the minimization subroutine SBMIN by Conn, Gould and Toint.

!  Nick Gould, for CGT productions.
!  October 1991.


! ---------------------------------------------------------------------




! ---------------------------------------------------------------------



! ---------------------------------------------------------------------


! ---------------------------------------------------------------------

!  integer variables from the GLOBAL common block.


!  integer variables from the LOCAL common block.


!  Integer variables from the PRFCTS common block.


!  local variables.

      INTEGER :: I, J, IG, IFSTAT, IGSTAT
      REAL ( KIND = wp ) :: FTT, ONE, ZERO
      PARAMETER ( ZERO = 0.0_wp, ONE = 1.0_wp )

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
      END IF

!  Compute the objective and constraint function values.

!                            block if there are no constraints. 
      F = ZERO
      IF ( data%numcon > 0 ) THEN
         DO 210 IG = 1, data%ng
            IF ( data%KNDOFC( IG ) == 0 ) THEN
               IF ( data%GXEQX( IG ) ) THEN
                  F = F + data%GSCALE( IG ) * data%FT( IG )
               ELSE
                  F = F + data%GSCALE( IG ) * data%GVALS( IG )
               END IF
            ELSE
               IF ( data%GXEQX( IG ) ) THEN
                  C( data%KNDOFC( IG ) ) = &
                      data%GSCALE( IG ) * data%FT( IG )
               ELSE
                  C( data%KNDOFC( IG ) ) = &
                      data%GSCALE( IG ) * data%GVALS( IG )
                  END IF
               END IF
  210    CONTINUE
      ELSE

!  There are no constraints, so we need not check data%KNDOFC( IG ).

         DO 220 IG = 1, data%ng
            IF ( data%GXEQX( IG ) ) THEN
               F = F + data%GSCALE( IG ) * data%FT( IG )
            ELSE
               F = F + data%GSCALE( IG ) * data%GVALS( IG )
            END IF
  220    CONTINUE
      END IF

!  Update the counters for the report tool.

      data%nc2of = data%nc2of + 1
      data%nc2cf = data%nc2cf + data%pnc
      RETURN

!  end of CFN.

      END
