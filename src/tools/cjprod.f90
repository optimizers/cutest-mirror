! ( Last modified on 10 Sepc 2004 at 16:55:38 )
!  Correction: 10/Sep/2004: undeclared integers variables declared
      SUBROUTINE CJPROD( data, N, M, GOTJ, JTRANS, X, V, LV, R, LR )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: N, M, LV, LR
      LOGICAL :: GOTJ, JTRANS
      REAL ( KIND = wp ) :: X( N ), V( LV ), R( LR )

!  Compute the matrix-vector product between the Jacobian matrix
!  of the constraints (JTRANS = .FALSE.), or its transpose 
! (JTRANS = .TRUE.) for the problem, and a given vector P. 
!  The result is placed in R. If GOTJ is .TRUE. the first derivatives 
!  are assumed to have already been computed. If the user is unsure, 
!  set GOTJ = .FALSE. the first time a product is required with the 
!  Jacobian evaluated at X. X is not used if GOTJ = .TRUE.

!  Nick Gould, for GOT/CUTEr productions.
!  June, 2003.


! ---------------------------------------------------------------------




! ---------------------------------------------------------------------



! ---------------------------------------------------------------------


! ---------------------------------------------------------------------





!  Integer variables from the PRFCTS common block.


!  Local variables

      INTEGER :: I, IG, J, ICON, K, IG1, II
      INTEGER :: L, IEL, NVAREL, NIN
      INTEGER :: IFSTAT, IGSTAT
      REAL ( KIND = wp ) :: ZERO, ONE, FTT, PROD, SCALEE
      PARAMETER ( ZERO = 0.0_wp, ONE = 1.0_wp )
      EXTERNAL :: RANGE 
      IF ( data%numcon == 0 ) RETURN

!  Check input data.

      IF ( ( JTRANS .AND. LV < M ) .OR.  &
 ( .NOT. JTRANS .AND. LV < N ) ) THEN
         IF ( IOUT > 0 ) WRITE( IOUT, 2010 )
         STOP
      END IF
      IF ( ( JTRANS .AND. LR < N ) .OR.  &
 ( .NOT. JTRANS .AND. LR < M ) ) THEN
         IF ( IOUT > 0 ) WRITE( IOUT, 2020 )
         STOP
      END IF

!  There are non-trivial group functions.

      IF ( .NOT. GOTJ ) THEN
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
            END IF
   70    CONTINUE

!  Evaluate the group derivative values.

         IF ( .NOT. data%altriv ) CALL GROUP ( data%GVALS( 1 ), data%ng, &
               data%FT( 1 ), data%GPVALU( 1 ), data%ng, &
               data%ITYPEG( 1 ), data%ISTGP( 1 ), &
               data%ICALCF( 1 ), &
               data%lcalcg, data%ng1, data%lcalcg, data%lcalcg, data%lgpvlu, &
               .TRUE., IGSTAT )
      END IF

!  Ensure that there is sufficient space.

      IF ( data%lwk2 < N ) THEN
         IF ( IOUT > 0 ) WRITE( IOUT, 2000 )
         STOP
      END IF

!  Form the product r = J(transpose) v

      IF ( JTRANS ) THEN

!  Initialize R

         DO 110 I = 1, N
            R( I ) = ZERO
  110    CONTINUE

!  Consider the IG-th group.

         DO 190 IG = 1, data%ng
            ICON = data%KNDOFC( IG )
            IF ( ICON > 0 ) THEN
               IG1 = IG + 1


!  Compute the product of v(i) with the (scaled) group derivative

               IF ( data%GXEQX( IG ) ) THEN
                  PROD = V( ICON ) * data%GSCALE( IG )
               ELSE
                  PROD = V( ICON ) * data%GSCALE( IG ) * data%GVALS( data%ng + IG )
               END IF

!  Loop over the group's nonlinear elements.

               DO 150 II = data%ISTADG( IG ),  &
                           data%ISTADG( IG1 ) - 1
                  IEL = data%IELING( II )
                  K = data%INTVAR( IEL )
                  L = data%ISTAEV( IEL )
                  NVAREL = data%ISTAEV( IEL + 1 ) - L
                  SCALEE = data%ESCALE( II ) * PROD
                  IF ( data%INTREP( IEL ) ) THEN

!  The IEL-th element has an internal representation.

                     NIN = data%INTVAR( IEL + 1 ) - K
                     CALL RANGE ( IEL, .TRUE., data%FUVALS( K ), &
                                  data%WRK( 1 ), NVAREL, NIN, &
                                  data%ITYPEE( IEL ), &
                                  NIN, NVAREL )
!DIR$ IVDEP
                     DO 130 I = 1, NVAREL
                        J = data%IELVAR( L )
                        R( J ) = R( J ) + SCALEE * data%WRK( I )
                        L = L + 1
  130                CONTINUE
                  ELSE

!  The IEL-th element has no internal representation.

!DIR$ IVDEP
                     DO 140 I = 1, NVAREL
                        J = data%IELVAR( L )
                        R( J ) = R( J ) + SCALEE * data%FUVALS( K )
                        K = K + 1
                        L = L + 1
  140                CONTINUE
                  END IF
  150          CONTINUE

!  Include the contribution from the linear element.

!DIR$ IVDEP
               DO 160 K = data%ISTADA( IG ), &
                          data%ISTADA( IG1 ) - 1
                  J = data%ICNA( K )
                  R( J ) = R( J ) + data%A( K ) * PROD
  160          CONTINUE
            END IF
  190    CONTINUE

!  Form the product r = J v

      ELSE

!  Consider the IG-th group.

         DO 290 IG = 1, data%ng
            ICON = data%KNDOFC( IG )
            IF ( ICON > 0 ) THEN
               IG1 = IG + 1
               PROD = ZERO

!  Compute the first derivative of the group.

!  Loop over the group's nonlinear elements.

               DO 250 II = data%ISTADG( IG ),  &
                           data%ISTADG( IG1 ) - 1
                  IEL = data%IELING( II )
                  K = data%INTVAR( IEL )
                  L = data%ISTAEV( IEL )
                  NVAREL = data%ISTAEV( IEL + 1 ) - L
                  SCALEE = data%ESCALE( II )
                  IF ( data%INTREP( IEL ) ) THEN

!  The IEL-th element has an internal representation.

                     NIN = data%INTVAR( IEL + 1 ) - K
                     CALL RANGE ( IEL, .TRUE., data%FUVALS( K ), &
                                  data%WRK( 1 ), NVAREL, NIN, &
                                  data%ITYPEE( IEL ), &
                                  NIN, NVAREL )
!DIR$ IVDEP
                     DO 230 I = 1, NVAREL
                        PROD = PROD + V( data%IELVAR( L ) ) * &
                                       SCALEE * data%WRK( I )
                        L = L + 1
  230                CONTINUE
                  ELSE

!  The IEL-th element has no internal representation.

!DIR$ IVDEP
                     DO 240 I = 1, NVAREL
                        PROD = PROD + V( data%IELVAR( L ) ) * &
                                   SCALEE * data%FUVALS( K )
                        K = K + 1
                        L = L + 1
  240                CONTINUE
                  END IF
  250          CONTINUE

!  Include the contribution from the linear element.

!DIR$ IVDEP
               DO 260 K = data%ISTADA( IG ), &
                          data%ISTADA( IG1 ) - 1
                  PROD = PROD + V( data%ICNA( K ) ) * data%A( K )
  260          CONTINUE

!  Multiply the product by the (scaled) group derivative

               IF ( data%GXEQX( IG ) ) THEN
                  R( ICON ) = PROD * data%GSCALE( IG )
               ELSE
                  R( ICON ) = PROD * data%GSCALE( IG ) * data%GVALS( data%ng + IG )
               END IF
            END IF
  290    CONTINUE
      END IF

!  Update the counters for the report tool.

      IF ( .NOT. GOTJ ) THEN
         data%nc2og = data%nc2og + 1
         data%nc2cg = data%nc2cg + data%pnc
      END IF
      RETURN

! Non-executable statements.

 2000 FORMAT( ' ** SUBROUTINE CJPROD: Increase the size of WK ' )
 2010 FORMAT( ' ** SUBROUTINE CJPROD: Increase the size of V ' )
 2020 FORMAT( ' ** SUBROUTINE CJPROD: Increase the size of R ' )

!  end of CJPROD.

      END
