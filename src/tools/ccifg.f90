! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CCIFG ( data, status, n, icon, X, ci, GCI, GRAD )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n, icon
      INTEGER, INTENT( OUT ) :: status
      REAL ( KIND = wp ) :: ci
      REAL ( KIND = wp ) :: X( n ), GCI( n )
      LOGICAL :: GRAD

!  Evaluate constraint function icon and possibly its gradient,
!  for constraints initially written in Standard Input Format (SIF).  
!  The constraint gradient is stored as a dense vector in array GCI;
!  that is, GCI( j ) is the partial derivative of constraint icon
!  with respect to variable j.
! (Subroutine CSCIFG performs the same calculations for a sparse
!  constraint gradient vector.)

!  Ingrid Bongartz
!  September 1994.

!  local variables.

      INTEGER :: i, j, iel, k, ig, ii, ig1, l, ll, neling
      INTEGER :: nin, nvarel, nelow, nelup, istrgv, iendgv
      INTEGER :: ifstat, igstat
      LOGICAL :: NONTRV
      EXTERNAL :: RANGE 
      REAL ( KIND = wp ) :: ftt, one, zero, gi, scalee
      PARAMETER ( zero = 0.0_wp, one = 1.0_wp )

!  Return if there are no constraints.

      IF ( data%numcon == 0 ) RETURN

!  Check input parameters.

      IF ( icon <= 0 ) THEN
         IF ( data%out > 0 ) WRITE( data%out, 2000 )
         status = 2 ; RETURN
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
         IF ( data%out > 0 ) WRITE( data%out, 2000 )
         status = 2 ; RETURN
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

!  evaluate the element functions

      CALL ELFUN( data%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,          &
                  data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,          &
                  data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,           &
                  data%lelvar, data%lntvar, data%lstadh, data%lstep,           &
                  data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,          &
                  1, ifstat )

!  Compute the group argument value FTT.
!  Consider only the group associated with constraint ICON.

      ftt = - data%B( ig )

!  Include the contribution from the linear element 
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
         data%GVALS( ig, 1 ) = data%FT( ig )
         data%GVALS( ig, 2 ) = one

!  otherwise, evaluate group ig

      ELSE
        CALL GROUP( data%GVALS, data%ng, data%FT, data%GPVALU, data%ng,        &
                    data%ITYPEG, data%ISTGP, data%ICALCF, data%ltypeg,         &
                    data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,         &
                    .FALSE., igstat )
      END IF

!  Compute the constraint function value.

      IF ( data%GXEQX( ig ) ) THEN
         ci = data%GSCALE( ig ) * data%FT( ig )
      ELSE
         ci = data%GSCALE( ig ) * data%GVALS( ig, 1 )

!  Update the constraint function evaluation counter

      data%nc2cf = data%nc2cf + 1
      END IF
      IF ( GRAD ) THEN

!  Update the constraint gradient evaluation counter

      data%nc2cg = data%nc2cg + 1

!  evaluate the element function derivatives

        CALL ELFUN( data%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,        &
                    data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,        &
                    data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,         &
                    data%lelvar, data%lntvar, data%lstadh, data%lstep,         &
                    data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,        &
                    2, ifstat )

!  evaluate the group derivative

         IF ( .NOT. data%GXEQX( ig ) )                                         &
          CALL GROUP( data%GVALS, data%ng, data%FT, data%GPVALU, data%ng,      &
                      data%ITYPEG, data%ISTGP, data%ICALCF, data%ltypeg,       &
                      data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,       &
                      .TRUE., igstat )

!  Compute the gradient.  Initialize the gradient vector as zero.

         DO 120 j = 1, n
            GCI( j ) = zero
  120    CONTINUE

!  Consider only group IG.

         ig1 = ig + 1
         istrgv = data%IWORK( data%lstagv + ig )
         iendgv = data%IWORK( data%lstagv + ig1 ) - 1
         NONTRV = .NOT. data%GXEQX( ig )

!  Compute the first derivative of the group.

         gi = data%GSCALE( ig )
         IF ( NONTRV ) gi = gi  * data%GVALS( ig, 2 )

!  The group has nonlinear elements.

         IF ( nelow <= nelup ) THEN
            data%WRK( data%IWORK( data%lsvgrp + istrgv :                       &
                                  data%lsvgrp + iendgv ) ) = zero

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
                  CALL RANGE ( iel, .TRUE., data%FUVALS( k ),                  &
                               data%WRK( n + 1 ), nvarel, nin,                 &
                               data%ITYPEE( iel ), nin, nvarel )
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

               IF ( ll <= n ) GCI( ll ) = gi * data%WRK( ll )
  190       CONTINUE

!  The group has only linear elements.

         ELSE

!  Allocate a gradient.

!DIR$ IVDEP
            DO 210 k = data%ISTADA( ig ), data%ISTADA( ig1 ) - 1
               ll = data%ICNA( k )

!  Include contributions from the first n variables only.

               IF ( ll <= n ) GCI( ll ) = gi * data%A( k )
  210       CONTINUE
         END IF
      END IF
      status = 0
      RETURN

! Non-executable statements.

 2000 FORMAT( ' ** SUBROUTINE CCIFG: invalid constraint index icon ' )

!  end of CCIFG.

      END
