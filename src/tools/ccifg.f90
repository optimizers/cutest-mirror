! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CCIFG ( data, N, ICON, X, CI, GCI, GRAD )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: N, ICON
      REAL ( KIND = wp ) :: CI
      REAL ( KIND = wp ) :: X( N ), GCI( N )
      LOGICAL :: GRAD

!  Evaluate constraint function ICON and possibly its gradient,
!  for constraints initially written in Standard Input Format (SIF).  
!  The constraint gradient is stored as a dense vector in array GCI;
!  that is, GCI( j ) is the partial derivative of constraint ICON
!  with respect to variable j.
! (Subroutine CSCIFG performs the same calculations for a sparse
!  constraint gradient vector.)

!  Ingrid Bongartz
!  September 1994.


! ---------------------------------------------------------------------




! ---------------------------------------------------------------------



! ---------------------------------------------------------------------


! ---------------------------------------------------------------------

!  integer variables from the GLOBAL common block.


!  integer variables from the LOCAL common block.


!  Integer variables from the PRFCTS common block.


!  local variables.

      INTEGER :: I, J, IEL, K, IG, II, IG1, L, LL, NELING
      INTEGER :: NIN, NVAREL, NELOW, NELUP, ISTRGV, IENDGV
      INTEGER :: IFSTAT, IGSTAT
      LOGICAL :: NONTRV
!D    EXTERNAL           DSETVL, DSETVI, RANGE 
      REAL ( KIND = wp ) :: FTT, ONE, ZERO, GI, SCALEE
      PARAMETER ( ZERO = 0.0_wp, ONE = 1.0_wp )

!  Return if there are no constraints.

      IF ( data%numcon == 0 ) RETURN

!  Check input parameters.

      IF ( ICON <= 0 ) THEN
         WRITE( IOUT, 2000 )
         STOP
      END IF

!  Find group index IG of constraint ICON.

      IG = 0
      DO 70 I = 1, data%ng
         IF ( data%KNDOFC( I ) == ICON ) THEN
           IG = I
           GOTO 80
         END IF
   70 CONTINUE
   80 IF ( IG == 0 ) THEN
         WRITE( IOUT, 2000 )
         STOP
      END IF

!  Determine nonlinear elements in group IG.
!  Record their indices in data%IWORK( ICALCF ).

      NELOW = data%ISTADG( IG )
      NELUP = data%ISTADG( IG + 1 ) - 1
      NELING = NELUP - NELOW + 1
      J = NELOW - 1
      DO 10 I = 1, NELING
         J = J + 1
         data%ICALCF( I ) = data%IELING( J )
   10 CONTINUE

!  Evaluate the element functions.

      CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), NELING, &
                   data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                   data%IELVAR( 1 ), data%INTVAR( 1 ), &
                   data%ISTADH( 1 ), data%ISTEP( 1 ), &
                   data%ICALCF( 1 ),  &
                   data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                   data%lntvar, data%lintre, LFUVAL, data%lvscal, data%lepvlu,  &
                   1, IFSTAT )

!  Compute the group argument value FTT.
!  Consider only the group associated with constraint ICON.

      FTT = - data%B( IG )

!  Include the contribution from the linear element 
!  only if the variable belongs to the first N variables.

      DO 30 I = data%ISTADA( IG ), data%ISTADA( IG + 1 ) - 1
         J = data%ICNA( I )
         IF ( J <= N ) &
            FTT = FTT + data%A( I ) * X( J )
   30 CONTINUE

!  Include the contributions from the nonlinear elements.

      DO 60 I = NELOW, NELUP
         FTT = FTT + &
            data%ESCALE( I ) * data%FUVALS( data%IELING( I ) )
   60 CONTINUE
      data%FT( IG ) = FTT

!  If IG is a trivial group, record the function value and derivative.

      IF ( data%GXEQX( IG ) ) THEN
         data%GVALS( IG ) = data%FT( IG )
         data%GVALS( data%ng + IG ) = ONE
      ELSE

!  Otherwise, evaluate group IG.

         CALL GROUP ( data%GVALS( 1 ), data%ng, data%FT( 1 ), &
                      data%GPVALU( 1 ), 1, &
                      data%ITYPEG( 1 ), data%ISTGP( 1 ), &
                      IG, data%lcalcg, data%ng1, data%lcalcg, data%lcalcg, data%lgpvlu, &
                      .FALSE., IGSTAT )
      END IF

!  Compute the constraint function value.

      IF ( data%GXEQX( IG ) ) THEN
         CI = data%GSCALE( IG ) * data%FT( IG )
      ELSE
         CI = data%GSCALE( IG ) * data%GVALS( IG )

!  Update the constraint function evaluation counter

      data%nc2cf = data%nc2cf + 1
      END IF
      IF ( GRAD ) THEN

!  Update the constraint gradient evaluation counter

      data%nc2cg = data%nc2cg + 1

!  Evaluate the element function derivatives.

         CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), NELING, &
                      data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                      data%IELVAR( 1 ), data%INTVAR( 1 ), &
                      data%ISTADH( 1 ), data%ISTEP( 1 ), &
                      data%ICALCF( 1 ),  &
                      data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                      data%lntvar, data%lintre, LFUVAL, data%lvscal, data%lepvlu,  &
                      2, IFSTAT )

!  Evaluate the group derivative.

         IF ( .NOT. data%GXEQX( IG ) ) &
            CALL GROUP ( data%GVALS( 1 ), data%ng, data%FT( 1 ), &
                         data%GPVALU( 1 ), 1, &
                         data%ITYPEG( 1 ), data%ISTGP( 1 ), &
                         IG, data%lcalcg, data%ng1, data%lcalcg, data%lcalcg, data%lgpvlu, &
                         .TRUE., IGSTAT )

!  Compute the gradient.  Initialize the gradient vector as zero.

         DO 120 J = 1, N
            GCI( J ) = ZERO
  120    CONTINUE

!  Consider only group IG.

         IG1 = IG + 1
         ISTRGV = data%IWORK( data%lstagv + IG )
         IENDGV = data%IWORK( data%lstagv + IG1 ) - 1
         NONTRV = .NOT. data%GXEQX( IG )

!  Compute the first derivative of the group.

         GI = data%GSCALE( IG )
         IF ( NONTRV ) GI = GI  * data%GVALS( data%ng + IG )

!  The group has nonlinear elements.

         IF ( NELOW <= NELUP ) THEN
      CALL DSETVI( IENDGV - ISTRGV + 1, data%WRK( 1 ), &
                         data%IWORK( data%lsvgrp + ISTRGV ), ZERO )

!  Loop over the group's nonlinear elements.

            DO 150 II = NELOW, NELUP
               IEL = data%IELING( II )
               K = data%INTVAR( IEL )
               L = data%ISTAEV( IEL )
               NVAREL = data%ISTAEV( IEL + 1 ) - L
               SCALEE = data%ESCALE( II )
               IF ( data%INTREP( IEL ) ) THEN

!  The IEL-th element has an internal representation.

                  NIN = data%INTVAR( IEL + 1 ) - K
                  CALL RANGE ( IEL, .TRUE., data%FUVALS( K ), &
                               data%WRK( N + 1 ), NVAREL, NIN, &
                               data%ITYPEE( IEL ), &
                               NIN, NVAREL )
!DIR$ IVDEP
                  DO 130 I = 1, NVAREL
                     J = data%IELVAR( L )
                     data%WRK( J ) = data%WRK( J ) + &
                                        SCALEE * data%WRK( N + I )
                     L = L + 1
  130             CONTINUE
               ELSE

!  The IEL-th element has no internal representation.

!DIR$ IVDEP
                  DO 140 I = 1, NVAREL
                     J = data%IELVAR( L )
                     data%WRK( J ) = data%WRK( J ) + &
                                        SCALEE * data%FUVALS( K )
                     K = K + 1
                     L = L + 1
  140             CONTINUE
               END IF
  150       CONTINUE

!  Include the contribution from the linear element.

!DIR$ IVDEP
            DO 160 K = data%ISTADA( IG ), &
                               data%ISTADA( IG1 ) - 1
               J = data%ICNA( K )
               data%WRK( J ) = data%WRK( J ) + data%A( K )
  160       CONTINUE

!  Allocate a gradient.

!DIR$ IVDEP
            DO 190 I = ISTRGV, IENDGV
               LL = data%IWORK( data%lsvgrp + I )

!  Include contributions from the first N variables only.

               IF ( LL <= N ) GCI( LL ) = GI * data%WRK( LL )
  190       CONTINUE

!  The group has only linear elements.

         ELSE

!  Allocate a gradient.

!DIR$ IVDEP
            DO 210 K = data%ISTADA( IG ), data%ISTADA( IG1 ) - 1
               LL = data%ICNA( K )

!  Include contributions from the first N variables only.

               IF ( LL <= N ) GCI( LL ) = GI * data%A( K )
  210       CONTINUE
         END IF
      END IF
      RETURN

! Non-executable statements.

 2000 FORMAT( ' ** SUBROUTINE CCIFG: invalid constraint index ICON ' )

!  end of CCIFG.

      END
