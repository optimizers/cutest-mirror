! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CGRDH ( data, n, m, X, GRLAGF, lv, V, G, &
                         JTRANS, lcjac1, lcjac2, CJAC, lh1, H )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n, m, lv, lh1, lcjac1, lcjac2
      LOGICAL :: GRLAGF, JTRANS
      REAL ( KIND = wp ) :: X ( n ), G ( n ), V ( lv ), &
                         CJAC ( lcjac1, lcjac2 ), &
                         H ( lh1, n )

!  Compute both the gradients of the objective, or Lagrangian, and
!  general constraint functions and the Hessian matrix of the
!  Lagrangian function of a problem initially written in Standard
!  Input Format (SIF).

!  G	 is an array which gives the value of the gradient of
!	 the objective function evaluated at X (GRLAGF = .FALSE.)
!        of of the Lagrangian function evaluated at X and V
! (GRLAGF = .TRUE.),

!  CJAC	 is a two-dimensional array of dimension ( lcjac1, lcjac2 )
!	 which gives the value of the Jacobian matrix of the
!	 constraint functions, or its transpose, evaluated at X.
!	 If JTRANS is .TRUE., the i,j-th component of the array
!        will contain the i-th derivative of the j-th constraint
!        function. Otherwise, if JTRANS is .FALSE., the i,j-th
!        component of the array will contain the j-th derivative
!        of the i-th constraint function.

!  H     is a two-dimensional array which gives the value of the
!        Hessian matrix of the Lagrangian function evaluated at
!        X and V. The i,j-th component of the array will contain the
!        derivative with respect to variables X(i) and X(j).

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


!  Integer variables from the PRFCTS common block.


!  Local variables

      INTEGER :: lh, lwkh, liwkh, lirnh, ljcnh, icon
      INTEGER :: lnxtrw, linxtr, inform, nnzh, iendgv
      INTEGER :: i, j, iel, k, ig, ii, ig1, l, jj, ll, lih
      INTEGER :: nin, nvarel, nelow, nelup, istrgv
      INTEGER :: ifstat, igstat
      LOGICAL :: NONTRV
!D    EXTERNAL           SETVL, SETVI, RANGE 
      REAL ( KIND = wp ) :: ftt, one, zero, gi, scalee, gii
      PARAMETER ( zero = 0.0_wp, one = 1.0_wp )

!  Check input parameters.

!                            dimension-checking.
      IF ( data%numcon > 0 ) THEN
         IF ( JTRANS ) THEN
            IF ( lcjac1 < n .OR. lcjac2 < m ) THEN
               IF ( lcjac1 < n ) WRITE( iout, 2020 )
               IF ( lcjac2 < m ) WRITE( iout, 2030 )
               STOP
            END IF
         ELSE
            IF ( lcjac1 < m .OR. lcjac2 < n ) THEN
               IF ( lcjac1 < m ) WRITE( iout, 2020 )
               IF ( lcjac2 < n ) WRITE( iout, 2030 )
               STOP
            END IF
         END IF
      END IF
      IF ( lh1 < n ) THEN
         WRITE( iout, 2040 )
         STOP
      END IF

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

!  evaluate the element function values.

      CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), data%nel, &
                   data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                   data%IELVAR( 1 ), data%INTVAR( 1 ), &
                   data%ISTADH( 1 ), data%ISTEP( 1 ), &
                   data%ICALCF( 1 ),  &
                   data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                   data%lntvar, data%lintre, lfuval, data%lvscal, data%lepvlu,  &
                   3, ifstat )

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

         IF ( data%GXEQX( ig ) ) THEN
            data%GVALS( data%ng + ig ) = one
            data%GVALS( 2 * data%ng + ig ) = zero
         END IF
   40 CONTINUE

!  evaluate the group derivative values.

      IF ( .NOT. data%altriv ) CALL GROUP ( data%GVALS( 1 ), data%ng, &
            data%FT( 1 ), data%GPVALU( 1 ), data%ng, &
            data%ITYPEG( 1 ), data%ISTGP( 1 ), &
            data%ICALCF( 1 ), &
            data%lcalcg, data%ng1, data%lcalcg, data%lcalcg, data%lgpvlu, &
            .TRUE., igstat )

!  Define the real work space needed for ELGRD.
!  Ensure that there is sufficient space.

      IF ( data%lwk2 < data%ng ) THEN
         IF ( iout > 0 ) WRITE( iout, 2000 )
         STOP
      END IF

!  For unconstrained problems, skip construction of data%WRK( ig ) 
!  and skip specialized construction of gradient and Jacobian.  
!  Call ELGRD instead.

      IF ( data%numcon > 0 ) THEN

!  Change the group weightings to include the contributions from
!  the Lagrange multipliers.

         DO 80 ig = 1, data%ng
            i = data%KNDOFC( ig )
            IF ( i == 0 ) THEN
               data%WRK( ig ) = data%GSCALE( ig )
            ELSE
               data%WRK( ig ) = data%GSCALE( ig ) * V( i )
            END IF
   80    CONTINUE

!  Compute the gradient values. Initialize the gradient and
!  Jacobian (or its transpose) as zero.

         DO 120 j = 1, n
            G( j ) = zero
            DO 110 i = 1, m
               IF ( JTRANS ) THEN
                  CJAC( j, i ) = zero
               ELSE
                  CJAC( i, j ) = zero
               END IF
  110       CONTINUE
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
            gii = data%WRK( ig )
            IF ( NONTRV ) THEN
               gi = gi  * data%GVALS( data%ng + ig )
               gii = gii * data%GVALS( data%ng + ig )
            END IF

!  This is the first gradient evaluation or the group has nonlinear
!  elements.

            IF ( data%firstg .OR. nelow <= nelup ) THEN
      CALL SETVI( iendgv - istrgv + 1, data%WRK( 1 ), &
                            data%IWORK( data%lsvgrp + istrgv ), zero )

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
                     G( ll ) = G( ll ) + gi * data%WRK( ll )

!  The group defines a constraint.

                  ELSE
                     IF ( JTRANS ) THEN
                        CJAC( ll, icon ) = gi * data%WRK( ll )
                     ELSE
                        CJAC( icon, ll ) = gi * data%WRK( ll )
                     END IF
                     IF ( GRLAGF ) &
                        G( ll ) = G( ll ) + gii * data%WRK( ll )
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

!  Allocate a gradient.

!DIR$ IVDEP
               DO 210 k = data%ISTADA( ig ), data%ISTADA( ig1 ) - 1
                  ll = data%ICNA( k )

!  The group belongs to the objective function.

                  IF ( icon == 0 ) THEN
                     G( ll ) = G( ll ) + gi * data%A( k )

!  The group defines a constraint.

                  ELSE
                     IF ( JTRANS ) THEN
                        CJAC( ll, icon ) =            gi * data%A( k )
                     ELSE
                        CJAC( icon, ll ) =            gi * data%A( k )
                     END IF
                     IF ( GRLAGF ) G( ll ) = G( ll ) + gii * data%A( k )
                  END IF
  210          CONTINUE

!  The group is non-trivial; increment the starting addresses for
!  the groups used by each variable in the (unchanged) linear
!  element to avoid resetting the nonzeros in the jacobian.

               IF ( NONTRV ) THEN
!DIR$ IVDEP
                  DO 220 i = istrgv, iendgv
                     ll = data%IWORK( data%lsvgrp + i )
                     data%IWORK( data%lstajc + ll ) = data%IWORK( data%lstajc + ll ) + 1
  220             CONTINUE
               END IF
            END IF
  290    CONTINUE

!  Reset the starting addresses for the lists of groups using
!  each variable to their values on entry.

         DO 300 i = n, 2, - 1
            data%IWORK( data%lstajc + i ) = data%IWORK( data%lstajc + i - 1 )
  300    CONTINUE
         data%IWORK( data%lstajc + 1 ) = 1
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

!  Store the gradient value.

         DO 400 i = 1, n
            G( i ) = data%FUVALS( data%lggfx + i )
  400    CONTINUE
      END IF
      data%firstg = .FALSE.

!  Define the real work space needed for ASMBL.
!  Ensure that there is sufficient space.

!                            for unconstrained problems.
      IF ( data%numcon > 0 ) THEN
         lwkh = data%lwk2 - n - 3 * data%maxsel - data%ng
      ELSE
         lwkh = data%lwk2 - n - 3 * data%maxsel
      END IF
      IF ( lwkh <= 0 ) THEN
         IF ( iout > 0 ) WRITE( iout, 2000 )
         STOP
      END IF

!  Define the integer work space needed for ASMBL.
!  Ensure that there is sufficient space.

      liwkh = data%liwk2 - n
      lh = MIN( lwkh, ( liwkh - 3 * n ) / 4 )
      linxtr = lh + n
      IF ( lh <= 0 ) THEN
         IF ( iout > 0 ) WRITE( iout, 2010 )
         STOP
      END IF

!  Set starting addresses for partitions of the integer workspace.

      lih = lh
      lirnh = 0
      ljcnh = lirnh + lih
      lnxtrw = ljcnh + lih
      DO 310 i = 1, n
         data%IVAR( i ) = i
  310 CONTINUE

!  Assemble the Hessian.

      IF ( data%numcon > 0 ) THEN
      CALL ASMBL( n, data%ng, data%maxsel, n, lh, lih, nnzh, &
                   n, data%IVAR( 1), data%ISTADH( 1 ), data%lstadh, &
                   data%ICNA( 1 ), data%licna, &
                   data%ISTADA( 1 ), data%lstada, data%INTVAR( 1 ), data%lntvar, &
                   data%IELVAR( 1 ), data%lelvar, data%IELING( 1 ), data%leling, &
                   data%ISTADG( 1 ), data%lstadg, data%ISTAEV( 1 ), data%lstaev, &
                   data%IWORK( data%lstagv + 1 ), data%lnstgv, data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, &
                   data%IWORK( data%lsend + lirnh + 1 ), data%IWORK( data%lsend + ljcnh + 1 ), &
                   data%IWORK( data%lsend + lnxtrw + 1 ), linxtr, &
                   data%IWORK( data%lsend + liwkh + 1 ), n, &
                   data%A( 1 ), data%la, data%FUVALS, data%lnguvl, data%FUVALS, data%lnhuvl, &
                   data%GVALS( data%ng + 1 ), data%GVALS( 2 * data%ng + 1 ), &
                   data%WRK( 1 ), data%ESCALE( 1 ), data%lescal, &
                   data%WRK( data%ng + 1 ), data%WRK( lwkh + data%ng + 1 ), &
                   data%lwk2 - lwkh, data%GXEQX( 1 ), &
                   data%lgxeqx, data%INTREP( 1 ), data%lintre, &
                   data%ITYPEE( 1 ), data%lintre, &
                   RANGE, 1, iout, .FALSE., i, inform, &
                   .FALSE., .FALSE. )
      ELSE
      CALL ASMBL( n, data%ng, data%maxsel, n, lh, lih, nnzh, &
                   n, data%IVAR( 1), data%ISTADH( 1 ), data%lstadh, &
                   data%ICNA( 1 ), data%licna, &
                   data%ISTADA( 1 ), data%lstada, data%INTVAR( 1 ), data%lntvar, &
                   data%IELVAR( 1 ), data%lelvar, data%IELING( 1 ), data%leling, &
                   data%ISTADG( 1 ), data%lstadg, data%ISTAEV( 1 ), data%lstaev, &
                   data%IWORK( data%lstagv + 1 ), data%lnstgv, data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, &
                   data%IWORK( data%lsend + lirnh + 1 ), data%IWORK( data%lsend + ljcnh + 1 ), &
                   data%IWORK( data%lsend + lnxtrw + 1 ), linxtr, &
                   data%IWORK( data%lsend + liwkh + 1 ), n, &
                   data%A( 1 ), data%la, data%FUVALS, data%lnguvl, data%FUVALS, data%lnhuvl, &
                   data%GVALS( data%ng + 1 ), data%GVALS( 2 * data%ng + 1 ), &
                   data%GSCALE( 1 ), data%ESCALE( 1 ), data%lescal, &
                   data%WRK( 1 ), data%WRK( lwkh + 1 ), &
                   data%lwk2 - lwkh, data%GXEQX( 1 ), &
                   data%lgxeqx, data%INTREP( 1 ), data%lintre, &
                   data%ITYPEE( 1 ), data%lintre, &
                   RANGE, 1, iout, .FALSE., i, inform, &
                   .FALSE., .FALSE. )
      END IF

!  Check that there is sufficient integer workspace.

      IF ( inform > 0 ) THEN
         IF ( iout > 0 ) WRITE( iout, 2010 )
         STOP
      END IF

!  Initialize the dense Hessian matrix.

      DO 420 j = 1, n
         DO 410 i = 1, n
            H( i, j ) = zero
  410    CONTINUE
  420 CONTINUE

!  Transfer the matrix from co-ordinate to dense storage and
!  symmetrize the martix.

      IF ( data%numcon > 0 ) THEN
         DO 430 k = 1, nnzh
            i = data%IWORK( data%lsend + lirnh + k )
            j = data%IWORK( data%lsend + ljcnh + k )
            H( i, j ) = data%WRK( data%ng + k )
            H( j, i ) = data%WRK( data%ng + k )
  430    CONTINUE
      ELSE
         DO 440 k = 1, nnzh
            i = data%IWORK( data%lsend + lirnh + k )
            j = data%IWORK( data%lsend + ljcnh + k )
            H( i, j ) = data%WRK( k )
            H( j, i ) = data%WRK( k )
  440    CONTINUE
      END IF

!  Update the counters for the report tool.

      data%nc2og = data%nc2og + 1
      data%nc2cg = data%nc2cg + data%pnc
      data%nc2oh = data%nc2oh + 1
      data%nc2ch = data%nc2ch + data%pnc
      RETURN

! Non-executable statements.

 2000 FORMAT( ' ** SUBROUTINE CGRDH: Increase the size of WK ' )
 2010 FORMAT( ' ** SUBROUTINE CGRDH: Increase the size of IWK ' )
 2020 FORMAT( ' ** SUBROUTINE CGRDH: Increase the leading dimension', &
              ' of CJAC ' )
 2030 FORMAT( ' ** SUBROUTINE CGRDH: Increase the second dimension', &
              ' of CJAC ' )
 2040 FORMAT( ' ** SUBROUTINE CGRDH: Increase the leading dimension', &
              ' of H ' )

!  end of CGRDH.

      END
