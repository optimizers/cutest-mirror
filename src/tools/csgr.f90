! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CSGR ( data, n, m, GRLAGF, lv, V, X, nnzj, &
                        lcjac, CJAC, INDVAR, INDFUN )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n, m, lv, nnzj, lcjac
      LOGICAL :: GRLAGF
      INTEGER :: INDVAR( LCJAC), INDFUN( lcjac )
      REAL ( KIND = wp ) :: X( n ), V ( lv ), CJAC ( lcjac )

!  Compute the gradients of the objective function and general
!  constraints of a function initially written in Standard
!  Input Format (SIF). The gradients are given in a sparse format.

!  CJAC	 is an array which gives the values of the nonzeros of the
!	 gradients of the objective, or Lagrangian, and general
!	 constraint functions evaluated  at X and V. The i-th entry
!	 of CJAC gives the value of the derivative with respect to
!	 variable INDVAR(i) of function INDFUN(i). INDFUN(i) = 0
!        indicates the objective function whenever GRLAGF is .FALSE.
!        or the Lagrangian function when GRLAGF is .TRUE., while
!        INDFUN(i) = j > 0 indicates the j-th general constraint
!        function.

!  Based on the minimization subroutine data%laNCELOT/SBMIN
!  by Conn, Gould and Toint.

!  Nick Gould, for CGT productions,
!  November 1991.


! ---------------------------------------------------------------------




! ---------------------------------------------------------------------



! ---------------------------------------------------------------------


! ---------------------------------------------------------------------

!  integer variables from the GLOBAL common block.


!  integer variables from the LOCAL common block.




!  local variables.

      INTEGER :: i, j, iel, k, ig, ii, ig1, l, jj, ll, icon
      INTEGER :: nin, nvarel, nelow, nelup, istrgv, iendgv
      INTEGER :: lgtemp, ifstat, igstat
      LOGICAL :: NONTRV
!D    EXTERNAL           SETVL, SETVI, RANGE 
      REAL ( KIND = wp ) :: ftt, one, zero, gi, scalee, gii
      PARAMETER ( zero = 0.0_wp, one = 1.0_wp )

!  there are non-trivial group functions.

      DO 10 i = 1, MAX( data%nelnum, data%ng )
        data%ICALCF( i ) = i
   10 CONTINUE

!  evaluate the element function values.

      CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), data%nelnum, &
                   data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                   data%IELVAR( 1 ), data%INTVAR( 1 ), &
                   data%ISTADH( 1 ), data%ISTEP( 1 ), &
                   data%ICALCF( 1 ),  &
                   data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                   data%lntvar, data%lintre, lfuval, data%lvscal, data%lepvlu,  &
                   1, ifstat )

!  evaluate the element function values.

      CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), data%nelnum, &
                   data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                   data%IELVAR( 1 ), data%INTVAR( 1 ), &
                   data%ISTADH( 1 ), data%ISTEP( 1 ), &
                   data%ICALCF( 1 ),  &
                   data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                   data%lntvar, data%lintre, lfuval, data%lvscal, data%lepvlu,  &
                   2, ifstat )

!  compute the group argument values ft.

      DO 40 ig = 1, data%ng
         ftt = - data%B( ig )

!  include the contribution from the linear element.

         DO 20 j = data%ISTADA( ig ), data%ISTADA( ig + 1 ) - 1
            ftt = ftt + data%A( j ) * X( data%ICNA( j ) )
   20    CONTINUE

!  include the contributions from the nonlinear elements.

         DO 30 j = data%ISTADG( ig ), data%ISTADG( ig + 1 ) - 1
            ftt = ftt + data%ESCALE( j ) * data%FUVALS( data%IELING( j ) )
   30    CONTINUE
         data%FT( ig ) = ftt

!  Record the derivatives of trivial groups.

         IF ( data%GXEQX( ig ) ) data%GVALS( data%ng + ig ) = one
   40 CONTINUE

!  evaluate the group derivative values.

      IF ( .NOT. data%altriv ) CALL GROUP ( data%GVALS( 1 ), data%ng, &
            data%FT( 1 ), data%GPVALU( 1 ), data%ng, &
            data%ITYPEG( 1 ), data%ISTGP( 1 ), &
            data%ICALCF( 1 ), &
            data%lcalcg, data%ng1, data%lcalcg, data%lcalcg, data%lgpvlu, &
            .TRUE., igstat )

!  Compute the gradient values. Initialize the gradient of the
!  objective function as zero.

      nnzj = 0
      IF ( data%numcon > 0 ) THEN
         lgtemp = WRK + n + data%maxsel
         DO 120 j = 1, n
            data%WRK( lgtemp + j ) = zero
  120    CONTINUE

!  Consider the IG-th group.

         DO 290 ig = 1, data%ng
            ig1 = ig + 1
            icon = data%KNDOFC( ig )
            istrgv = data%IWORK( data%lstagv + ig )
            iendgv = data%IWORK( data%lstagv + ig1 ) - 1
            nelow = data%ISTADG( ig )
            nelup = data%ISTADG( ig1 ) - 1
            NONTRV = .NOT. data%GXEQX( ig )

!  Compute the first derivative of the group.

            gi = data%GSCALE( ig )
            IF ( icon == 0 ) THEN
               gii = gi
            ELSE
               IF ( GRLAGF ) gii = gi * V( data%KNDOFC( ig ) )
            END IF
            IF ( NONTRV ) THEN
               gi = gi  * data%GVALS( data%ng + ig )
               IF ( GRLAGF ) gii = gii * data%GVALS( data%ng + ig )
            END IF
      CALL SETVI( iendgv - istrgv + 1, data%WRK( 1 ), &
                         data%IWORK( data%lsvgrp + istrgv ), zero )

!  This is the first gradient evaluation or the group has nonlinear
!  elements.

            IF ( data%firstg .OR. nelow <= nelup ) THEN

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

!  The group belongs to the objective function.

                  IF ( icon == 0 ) THEN
                     data%WRK( lgtemp + ll ) = data%WRK( lgtemp + ll ) + &
                                         gi * data%WRK( ll )

!  The group defines a constraint.

                  ELSE
                     nnzj = nnzj + 1
                     IF ( nnzj <= lcjac ) THEN
                        CJAC ( nnzj ) = gi * data%WRK( ll )
                        INDFUN( nnzj ) = icon
                        INDVAR( nnzj ) = ll
                     END IF
                     IF ( GRLAGF ) &
                        data%WRK( lgtemp + ll ) = data%WRK( lgtemp + ll ) + &
                                            gii * data%WRK( ll )
                  END IF

!  If the group is non-trivial, also store the nonzero entries of the
!  gradient of the function in GRJAC.

                  IF ( NONTRV ) THEN
                     jj = data%IWORK( data%lstajc + ll )
                     data%FUVALS( data%lgrjac + jj ) = data%WRK( ll )

!  Increment the address for the next nonzero in the column of
!  the jacobian for variable LL.

                     data%IWORK( data%lstajc + ll ) = jj + 1
                  END IF
  190          CONTINUE

!  This is not the first gradient evaluation and there is only a linear
!  element.

            ELSE
!                             linear element improved. 43 lines replace 40

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

!  The group belongs to the objective function.

                  IF ( icon == 0 ) THEN
                     data%WRK( lgtemp + ll ) = data%WRK( lgtemp + ll ) + &
                                         gi * data%WRK( ll )

!  The group defines a constraint.

                  ELSE
                     nnzj = nnzj + 1
                     IF ( nnzj <= lcjac ) THEN
                        CJAC ( nnzj ) = gi * data%WRK( ll )
                        INDFUN( nnzj ) = icon
                        INDVAR( nnzj ) = ll
                     END IF
                     IF ( GRLAGF ) &
                        data%WRK( lgtemp + ll ) = data%WRK( lgtemp + ll ) + &
                                            gii * data%WRK( ll )
                  END IF

!  Increment the address for the next nonzero in the column of
!  the jacobian for variable LL.

                  IF ( NONTRV ) THEN
                     jj = data%IWORK( data%lstajc + ll )
                     data%IWORK( data%lstajc + ll ) = jj + 1
                  END IF
  220          CONTINUE
            END IF
  290    CONTINUE

!  Reset the starting addresses for the lists of groups using
!  each variable to their values on entry.

         DO 300 i = n, 2, - 1
            data%IWORK( data%lstajc + i ) = data%IWORK( data%lstajc + i - 1 )
  300    CONTINUE
         data%IWORK( data%lstajc + 1 ) = 1

!  Transfer the gradient of the objective function to the sparse
!  storage scheme.

         DO 310 i = 1, n
!           IF ( data%WRK( lgtemp + i ) /= zero ) THEN
               nnzj = nnzj + 1
               IF ( nnzj <= lcjac ) THEN
                  CJAC ( nnzj ) = data%WRK( lgtemp + i )
                  INDFUN( nnzj ) = 0
                  INDVAR( nnzj ) = i
               END IF
!           END IF
  310    CONTINUE
      ELSE

!  Compute the gradient value.

      CALL ELGRD( n, data%ng, data%firstg, data%ICNA( 1 ), data%licna, &
                      data%ISTADA( 1 ), data%lstada, data%IELING( 1 ), &
                      data%leling, data%ISTADG( 1 ), data%lstadg, &
                      data%ITYPEE( 1 ), data%lintre, &
                      data%ISTAEV( 1 ), data%lstaev, data%IELVAR( 1 ), &
                      data%lelvar, data%INTVAR( 1 ), data%lntvar, &
                      data%IWORK( data%lsvgrp + 1 ), &
                      data%lnvgrp, data%IWORK( data%lstajc + 1 ), data%lnstjc, &
                      data%IWORK( data%lstagv + 1 ), data%lnstgv, data%A( 1 ), data%la, &
                      data%GVALS( data%ng + 1 ), data%lgvals, &
                      data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ), &
                      data%GSCALE( 1 ), data%lgscal, &
                      data%ESCALE( 1 ), data%lescal, data%FUVALS( data%lgrjac + 1 ), &
                      data%lngrjc, data%WRK( 1 ), data%WRK( n + 1 ), data%maxsel, &
                      data%GXEQX( 1 ), data%lgxeqx, &
                      data%INTREP( 1 ), data%lintre, RANGE )

!  Transfer the gradient of the objective function to the sparse
!  storage scheme.

         DO 400 i = 1, n
!           IF ( data%FUVALS( data%lggfx + i ) /= zero ) THEN
               nnzj = nnzj + 1
               IF ( nnzj <= lcjac ) THEN
                  CJAC ( nnzj ) = data%FUVALS( data%lggfx + i )
                  INDFUN( nnzj ) = 0
                  INDVAR( nnzj ) = i
               END IF
!           END IF
  400    CONTINUE
      END IF
      data%firstg = .FALSE.

!  Verify that the Jacobian can fit in the alloted space

      IF ( nnzj > lcjac ) THEN
         IF ( iout > 0 ) WRITE( iout, 2000 ) nnzj - lcjac 
         STOP
      END IF

!  Update the counters for the report tool.

      data%nc2og = data%nc2og + 1
      data%nc2cg = data%nc2cg + data%pnc
      RETURN

!  Non-executable statements.

 2000 FORMAT( /, ' ** SUBROUTINE CSGR: array length lcjac too small.',  &
              /, ' -- Minimization abandoned.', &
              /, ' -- Increase the parameter lcjac by at least ', I8, &
                 ' and restart.' )

!  end of CSGR.

      END
