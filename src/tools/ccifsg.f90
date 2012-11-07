! ( Last modified on 10 Sepc 2004 at 17:01:38 )
!  Correction: 10/Sep/2004: Incorrect argument to RANGES removed
      SUBROUTINE CCIFSG ( data, n, icon, X, ci, nnzgci, lgci, GCI, &
                          INDVAR, GRAD )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n, icon, nnzgci, lgci
      LOGICAL :: GRAD
      INTEGER :: INDVAR( lgci )
      REAL ( KIND = wp ) :: ci
      REAL ( KIND = wp ) :: X( n ), GCI( lgci )

!  Evaluate constraint function icon and possibly its gradient,
!  for constraints initially written in Standard Input Format (SIF).
!  The constraint gradient is stored as a sparse vector in array GCI.
!  The j-th entry of GCI gives the value of the partial derivative
!  of constraint icon with respect to variable INDVAR( j ).
!  The number of nonzeros in vector GCI is given by NNZGCI.
! (Subroutine CCIFG performs the same calculations for a dense
!  constraint gradient vector.)

!  Ingrid Bongartz
!  September 1994.


! ---------------------------------------------------------------------




! ---------------------------------------------------------------------


! ---------------------------------------------------------------------


! ---------------------------------------------------------------------


!  integer variables from the GLOBAL common block.


!  integer variables from the LOCAL common block.


!  integer variables from the PRFCTS common block.


!  local variables.

      INTEGER :: i, j, iel, k, ig, ii, ig1, l, ll, neling
      INTEGER :: nin, nvarel, nelow, nelup, istrgv, iendgv
      INTEGER :: ifstat, igstat
      LOGICAL :: NONTRV
!D    EXTERNAL           SETVL, SETVI, RANGE 
      REAL ( KIND = wp ) :: ftt, one, zero, gi, scalee
      PARAMETER ( zero = 0.0_wp, one = 1.0_wp )

!  Return if there are no constraints.

      IF ( data%numcon == 0 ) RETURN

!  Check input parameters.

      IF ( icon <= 0 ) THEN
         WRITE( iout, 2000 )
         STOP
      END IF

!  Find group index ig of constraint ICON.

      ig = 0
      DO 70 i = 1, data%ng
         IF ( data%KNDOFC( i ) == icon ) THEN
           ig = i
           GOTO 80
         END IF
   70 CONTINUE
   80 IF ( ig == 0 ) THEN
         WRITE( iout, 2000 )
         STOP
      END IF

!  Determine nonlinear elements in group IG.
!  Record their indices in data%IWORK( ICALCF ).

      nelow = data%ISTADG( ig )
      nelup = data%ISTADG( ig + 1 ) - 1
      neling = nelup - nelow + 1
      j = nelow - 1
      DO 10 i = 1, neling
         j = j + 1
         data%ICALCF( i ) = data%IELING( j )
   10 CONTINUE

!  Evaluate the element function values.

      CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), neling, &
                   data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                   data%IELVAR( 1 ), data%INTVAR( 1 ), &
                   data%ISTADH( 1 ), data%ISTEP( 1 ), &
                   data%ICALCF( 1 ),  &
                   data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                   data%lntvar, data%lintre, lfuval, data%lvscal, data%lepvlu,  &
                   1, ifstat )

!  Compute the group argument value FTT.
!  Consider only the group associated with constraint ICON.

      ftt = - data%B( ig )

!  Include contributions from the linear element 
!  only if the variable belongs to the first n variables.

      DO 30 i = data%ISTADA( ig ), data%ISTADA( ig + 1 ) - 1
         j = data%ICNA( i )
         IF ( j <= n ) &
            ftt = ftt + data%A( i ) * X( j )
   30 CONTINUE

!  Include the contributions from the nonlinear elements.

      DO 60 i = nelow, nelup
         ftt = ftt + &
                data%ESCALE( i ) * data%FUVALS( data%IELING( i ) )
   60 CONTINUE
      data%FT( ig ) = ftt

!  If ig is a trivial group, record the function value and derivative.

      IF ( data%GXEQX( ig ) ) THEN
         data%GVALS( ig ) = data%FT( ig )
         data%GVALS( data%ng + ig ) = one
      ELSE

!  Otherwise, evaluate group IG.

         CALL GROUP ( data%GVALS( 1 ), data%ng, data%FT( 1 ), &
                      data%GPVALU( 1 ), 1, &
                      data%ITYPEG( 1 ), data%ISTGP( 1 ), &
                      ig, data%lcalcg, data%ng1, data%lcalcg, data%lcalcg, data%lgpvlu, &
                      .FALSE., igstat )
      END IF

!  Compute the constraint function value.

      IF ( data%GXEQX( ig ) ) THEN
         ci = data%GSCALE( ig ) * data%FT( ig )
      ELSE 
         ci = data%GSCALE( ig ) * data%GVALS( ig )
      END IF
      IF ( GRAD ) THEN

!  Evaluate the element function derivatives.

         CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), neling, &
                      data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                      data%IELVAR( 1 ), data%INTVAR( 1 ), &
                      data%ISTADH( 1 ), data%ISTEP( 1 ), &
                      data%ICALCF( 1 ),  &
                      data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                      data%lntvar, data%lintre, lfuval, data%lvscal, data%lepvlu,  &
                      2, ifstat )

!  Evaluate the group derivative values.

         IF ( .NOT. data%GXEQX( ig ) ) &
            CALL GROUP ( data%GVALS( 1 ), data%ng, data%FT( 1 ), &
                         data%GPVALU( 1 ), 1, &
                         data%ITYPEG( 1 ), data%ISTGP( 1 ), &
                         ig, data%lcalcg, data%ng1, data%lcalcg, data%lcalcg, data%lgpvlu, &
                         .TRUE., igstat )

!  Compute the gradient values.  Initialize the gradient vector as zero.

         nnzgci = 0 
         DO 120 j = 1, lgci
            GCI( j ) = zero
  120    CONTINUE

!  Consider only group IG.

         ig1 = ig + 1
         istrgv = data%IWORK( data%lstagv + ig )
         iendgv = data%IWORK( data%lstagv + ig1 ) - 1
         NONTRV = .NOT. data%GXEQX( ig )

!  Compute the first derivative of the group.

         gi = data%GSCALE( ig )
         IF ( NONTRV ) gi = gi  * data%GVALS( data%ng + ig )
      CALL SETVI( iendgv - istrgv + 1, data%WRK( 1 ), &
                      data%IWORK( data%lsvgrp + istrgv ), zero )

!  The group has nonlinear elements.

         IF ( nelow <= nelup ) THEN

!  Loop over the group's nonlinear elements.

            DO 150 ii = nelow, nelup
               iel = data%IELING( ii )
               k = data%INTVAR( iel )
               l = data%ISTAEV( iel )
               nvarel = data%ISTAEV( iel + 1 ) - l
               scalee = data%ESCALE( ii )
               IF ( data%INTREP( iel ) ) THEN

!  The IEL-th element has an internal representation.

                  nin = data%INTVAR( iel + 1 ) - k
                  CALL RANGE ( iel, .TRUE., data%FUVALS( k ), &
                               data%WRK( n + 1 ), nvarel, nin, &
                               data%ITYPEE( iel ), &
                               nin, nvarel )
!DIR$ IVDEP
                  DO 130 i = 1, nvarel
                     j = data%IELVAR( l )
                     data%WRK( j ) = data%WRK( j ) + &
                                        scalee * data%WRK( n + i )
                     l = l + 1
  130             CONTINUE
               ELSE

!  The IEL-th element has no internal representation.

!DIR$ IVDEP
                  DO 140 i = 1, nvarel
                     j = data%IELVAR( l )
                     data%WRK( j ) = data%WRK( j ) + &
                                        scalee * data%FUVALS( k )
                     k = k + 1
                     l = l + 1
  140             CONTINUE
               END IF
  150       CONTINUE

!  Include the contribution from the linear element.

!DIR$ IVDEP
            DO 160 k = data%ISTADA( ig ), &
                               data%ISTADA( ig1 ) - 1
               j = data%ICNA( k )
               data%WRK( j ) = data%WRK( j ) + data%A( k )
  160       CONTINUE

!  Allocate a gradient.

!DIR$ IVDEP
            DO 190 i = istrgv, iendgv
               ll = data%IWORK( data%lsvgrp + i )

!  Include contributions from the first n variables only.

               IF ( ll <= n ) THEN
                  nnzgci = nnzgci + 1
                  GCI ( nnzgci ) = gi * data%WRK( ll )
                  INDVAR( nnzgci ) = ll
               END IF
  190       CONTINUE

!  The group has only linear elements.

         ELSE

!  Include the contribution from the linear element.

!DIR$ IVDEP
            DO 210 k = data%ISTADA( ig ),data%ISTADA( ig1 ) - 1
               j = data%ICNA( k )
               data%WRK( j ) = data%WRK( j ) + data%A( k )
  210       CONTINUE

!  Allocate a gradient.

!DIR$ IVDEP
            DO 220 i = istrgv, iendgv
               ll = data%IWORK( data%lsvgrp + i )

!  Include contributions from the first n variables only.

               IF ( ll <= n ) THEN
                  nnzgci = nnzgci + 1
                  GCI ( nnzgci ) = gi * data%WRK( ll )
                  INDVAR( nnzgci ) = ll
               END IF
  220       CONTINUE
         END IF
      END IF

!  Update the counters for the report tool.

      data%nc2cf = data%nc2cf + 1
      IF ( GRAD ) data%nc2cg = data%nc2cg + 1
      RETURN

! Non-executable statements.

 2000 FORMAT( ' ** SUBROUTINE CCIFSG: invalid constraint index icon ' )

!  end of CCIFSG.

      END
