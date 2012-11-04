! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CCFSG ( data, n, m, X, lc, C, nnzj, lcjac, CJAC, &
                         INDVAR, INDFUN, GRAD )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n, m, lc, nnzj, lcjac
      LOGICAL :: GRAD
      INTEGER :: INDVAR( lcjac ), INDFUN( lcjac )
      REAL ( KIND = wp ) :: X( n ), C( lc ), CJAC ( lcjac )

!  Compute the values of the constraint functions and their gradients
!  for constraints initially written in Standard Input Format (SIF).
!  The Jacobian must be stored in a sparse format.
! (Subroutine CCFG performs the same calculations for a dense Jacobian.)

!  CJAC  is an array which gives the values of the nonzeros of the
!        general constraint functions evaluated at X and V.
!        The i-th entry of CJAC gives the value of the derivative
!        with respect to variable INDVAR(i) of constraint function 
!        INDFUN(i) (i.e., INDFUN(i) = j > 0 indicates the j-th
!        general constraint function).

!  Based on the subroutines cfn.f and csgr.f by Nick Gould, which are
!  in turn based on the subroutine SBMIN by Conn, Gould and Toint.

!  Ingrid Bongartz
!  April 1992.


! ---------------------------------------------------------------------




! ---------------------------------------------------------------------



! ---------------------------------------------------------------------


! ---------------------------------------------------------------------

!  integer variables from the GLOBAL common block.


!  integer variables from the LOCAL common block.


!  integer variables from the PRFCTS common block.


!  local variables.

      INTEGER :: i, j, iel, k, ig, ii, ig1, l, ll, icon, icnt
      INTEGER :: nin, nvarel, nelow, nelup, istrgv, iendgv
      INTEGER :: llo, llwrk, ifstat, igstat
      LOGICAL :: NONTRV
!D    EXTERNAL           SETVL, SETVI, RANGE 
      REAL ( KIND = wp ) :: ftt, one, zero, gi, scalee
      PARAMETER ( zero = 0.0_wp, one = 1.0_wp )

      IF ( data%numcon == 0 ) RETURN

!  Must identify which elements are included in constraints.
!  Use logical work vector to keep track of elements already included.
!  First ensure there is sufficient room in LOGI.

      llo = gxeqx + data%ngng
      llwrk = llogic - llo
      IF ( llwrk < data%nelnum ) THEN
          IF ( iout > 0 ) WRITE( iout, 2010 ) data%nelnum - llwrk 
          STOP
      END IF
      DO 410 i = 1, data%nelnum
         data%LOGIC( i ) = .FALSE.
  410 CONTINUE

!  Now identify elements in first m constraint groups.

      icnt = 0
      DO 10 ig = 1, data%ng
         icon = data%KNDOFC( ig )
         IF ( icon > 0 .AND. icon <= m ) THEN
            nelow = data%ISTADG( ig )
            nelup = data%ISTADG( ig + 1 ) - 1
            DO 20 ii = nelow, nelup
               iel = data%IELING( ii )
               IF ( .NOT. data%LOGIC( iel ) ) THEN
                  data%LOGIC( iel ) = .TRUE.
                  icnt = icnt + 1
                  data%ICALCF( icnt ) = iel
               END IF
   20       CONTINUE
         END IF
   10 CONTINUE

!  Evaluate the element function values.

      CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), icnt, &
                   data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                   data%IELVAR( 1 ), data%INTVAR( 1 ), &
                   data%ISTADH( 1 ), data%ISTEP( 1 ), &
                   data%ICALCF( 1 ),  &
                   data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                   data%lntvar, data%lintre, lfuval, data%lvscal, data%lepvlu,  &
                   1, ifstat )
      IF ( GRAD ) THEN

!  Evaluate the element function derivatives.

         CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), icnt, &
                      data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                      data%IELVAR( 1 ), data%INTVAR( 1 ), &
                      data%ISTADH( 1 ), data%ISTEP( 1 ), &
                      data%ICALCF( 1 ),  &
                      data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                      data%lntvar, data%lintre, lfuval, data%lvscal, data%lepvlu,  &
                      2, ifstat )
      END IF

!  Compute the group argument values ft.

      DO 100 ig = 1, data%ng
         ftt = zero

!  Consider only those groups in the constraints.

         icon = data%KNDOFC( ig )
         IF ( icon > 0 .AND. icon <= m ) THEN
            ftt = - data%B( ig )

!  Include contributions from the linear element 
!  only if the variable belongs to the first n variables.

            DO 30 i = data%ISTADA( ig ), data%ISTADA( ig + 1 ) - 1
               j = data%ICNA( i )
               IF ( j <= n ) &
                  ftt = ftt + data%A( i ) * X( j )
   30       CONTINUE

!  Include the contributions from the nonlinear elements.

            DO 60 i = data%ISTADG( ig ), data%ISTADG( ig + 1 ) - 1
               ftt = ftt + &
                      data%ESCALE( i ) * data%FUVALS( data%IELING( i ) )
   60       CONTINUE

!  Record derivatives of trivial groups.

            IF ( data%GXEQX( ig ) ) data%GVALS( data%ng + ig ) = one 
         END IF
         data%FT( ig ) = ftt
  100 CONTINUE

!  Compute the group function values.

!  All group functions are trivial.

      IF ( data%altriv ) THEN
      CALL DCOPY( data%ng, data%FT( 1 ), 1, data%GVALS( 1 ), 1 )
      CALL SETVL( data%ng, data%GVALS( data%ng + 1 ), 1, one )
      ELSE

!  Evaluate the group function values.
!  Evaluate groups belonging to the first m constraints only.

         icnt = 0
         DO 400 ig = 1, data%ng
            icon = data%KNDOFC( ig )
            IF ( icon > 0 .AND. icon <= m ) THEN
               icnt = icnt + 1
               data%ICALCF( icnt ) = ig
            END IF 
  400    CONTINUE
         CALL GROUP ( data%GVALS( 1 ), data%ng, data%FT( 1 ), &
                      data%GPVALU( 1 ), icnt, &
                      data%ITYPEG( 1 ), data%ISTGP( 1 ), &
                      data%ICALCF( 1 ), &
                      data%lcalcg, data%ng1, data%lcalcg, data%lcalcg, data%lgpvlu, &
                      .FALSE., igstat )
      END IF

!  Compute the constraint function values.

      DO 110 ig = 1, data%ng
         i = data%KNDOFC( ig )
         IF ( i > 0 .AND. i <= m ) THEN
            IF ( data%GXEQX( ig ) ) THEN
               C( i ) = data%GSCALE( ig ) * data%FT( ig )
            ELSE 
               C( i ) = data%GSCALE( ig ) * data%GVALS( ig )
            END IF
         END IF
  110 CONTINUE
      IF ( GRAD ) THEN

!  Evaluate the group derivative values.

         IF ( .NOT. data%altriv ) CALL GROUP ( data%GVALS( 1 ), data%ng, &
               data%FT( 1 ), data%GPVALU( 1 ), icnt, &
               data%ITYPEG( 1 ), data%ISTGP( 1 ), &
               data%ICALCF( 1 ), &
               data%lcalcg, data%ng1, data%lcalcg, data%lcalcg, data%lgpvlu, &
               .TRUE., igstat )

!  Compute the gradient values.  Initialize the Jacobian as zero.

         nnzj = 0 
         DO 120 j = 1, lcjac
            CJAC( j ) = zero
  120    CONTINUE

!  Consider the IG-th group.

         DO 290 ig = 1, data%ng
            icon = data%KNDOFC( ig )

!  Consider only those groups in the first m constraints.

            IF ( icon == 0 .OR. icon > m ) GO TO 290
            ig1 = ig + 1
            istrgv = data%IWORK( data%lstagv + ig )
            iendgv = data%IWORK( data%lstagv + ig1 ) - 1
            nelow = data%ISTADG( ig )
            nelup = data%ISTADG( ig1 ) - 1
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
  130                CONTINUE
                  ELSE

!  The IEL-th element has no internal representation.

!DIR$ IVDEP
                     DO 140 i = 1, nvarel
                        j = data%IELVAR( l )
                        data%WRK( j ) = data%WRK( j ) + &
                                           scalee * data%FUVALS( k )
                        k = k + 1
                        l = l + 1
  140                CONTINUE
                  END IF
  150          CONTINUE

!  Include the contribution from the linear element.

!DIR$ IVDEP
               DO 160 k = data%ISTADA( ig ), &
                                  data%ISTADA( ig1 ) - 1
                  j = data%ICNA( k )
                  data%WRK( j ) = data%WRK( j ) + data%A( k )
  160          CONTINUE

!  Allocate a gradient.

!DIR$ IVDEP
               DO 190 i = istrgv, iendgv
                  ll = data%IWORK( data%lsvgrp + i )

!  Include contributions from the first n variables only.

                  IF ( ll <= n ) THEN
                     nnzj = nnzj + 1
                     IF ( nnzj <= lcjac ) THEN
                        CJAC ( nnzj ) = gi * data%WRK( ll )
                        INDFUN( nnzj ) = icon
                        INDVAR( nnzj ) = ll
                     END IF
                  END IF
  190          CONTINUE

!  The group has only linear elements.

            ELSE
!                             linear element improved. 26 lines replace 19

!  Include the contribution from the linear element.

!DIR$ IVDEP
               DO 210 k = data%ISTADA( ig ),data%ISTADA( ig1 ) - 1
                  j = data%ICNA( k )
                  data%WRK( j ) = data%WRK( j ) + data%A( k )
  210          CONTINUE

!  Allocate a gradient.

!DIR$ IVDEP
               DO 220 i = istrgv, iendgv
                  ll = data%IWORK( data%lsvgrp + i )

!  Include contributions from the first n variables only.

                  IF ( ll <= n ) THEN
                     nnzj = nnzj + 1
                     IF ( nnzj <= lcjac ) THEN
                        CJAC ( nnzj ) = gi * data%WRK( ll )
                        INDFUN( nnzj ) = icon
                        INDVAR( nnzj ) = ll
                     END IF
                  END IF
  220          CONTINUE
            END IF
  290    CONTINUE

!  Verify that the Jacobian can fit in the allotted space

         IF ( nnzj > lcjac ) THEN
            IF ( iout > 0 ) WRITE( iout, 2000 ) nnzj - lcjac 
            STOP
         END IF
      END IF

!  Update the counters for the report tool.

      data%nc2cf = data%nc2cf + data%pnc
      IF ( GRAD ) data%nc2cg = data%nc2cg + data%pnc
      RETURN

!  Non-executable statements.

 2000 FORMAT( /  ' ** SUBROUTINE CCFSG: array length lcjac too small.' &
              /  ' -- Minimization abandoned.', &
              /  ' -- Increase the parameter lcjac by at least ', I8, &
                 ' and restart.' )
 2010 FORMAT( /  ' ** SUBROUTINE CCFSG: array length llogic too small.' &
              /  ' -- Minimization abandoned.', &
              /  ' -- Increase the parameter llogic by at least ', I8, &
                 ' and restart.' )

!  end of CSCFG.

      END

