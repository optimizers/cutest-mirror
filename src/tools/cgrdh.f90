! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CGRDH( data, status, n, m, X, grlagf, Y, G,                   &
                        jtrans, lcjac1, lcjac2, CJAC, lh1, H )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n, m, lh1, lcjac1, lcjac2
      INTEGER, INTENT( OUT ) :: status
      LOGICAL :: grlagf, jtrans
      REAL ( KIND = wp ) :: X( n ), G( n ), Y( m ),                            &
                         CJAC( lcjac1, lcjac2 ), H( lh1, n )

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

!  Local variables

      INTEGER :: lh, lwkh, liwkh, lirnh, ljcnh, icon
      INTEGER :: lnxtrw, linxtr, inform, nnzh, iendgv
      INTEGER :: i, j, iel, k, ig, ii, ig1, l, jj, ll, lih
      INTEGER :: nin, nvarel, nelow, nelup, istrgv
      INTEGER :: ifstat, igstat
      LOGICAL :: NONTRV
      EXTERNAL :: RANGE 
      REAL ( KIND = wp ) :: ftt, gi, scalee, gii

!  check input parameters

      IF ( data%numcon > 0 ) THEN
        IF ( JTRANS ) THEN
          IF ( lcjac1 < n .OR. lcjac2 < m ) THEN
            IF ( lcjac1 < n .AND. data%out > 0 ) WRITE( data%out, 2020 )
            IF ( lcjac2 < m .AND. data%out > 0 ) WRITE( data%out, 2030 )
            status = 2 ; RETURN
          END IF
        ELSE
          IF ( lcjac1 < m .OR. lcjac2 < n ) THEN
            IF ( lcjac1 < m .AND. data%out > 0 ) WRITE( data%out, 2020 )
            IF ( lcjac2 < n .AND. data%out > 0 ) WRITE( data%out, 2030 )
            status = 2 ; RETURN
          END IF
        END IF
      END IF
      IF ( lh1 < n ) THEN
        IF ( data%out > 0 ) WRITE( data%out, 2040 )
        status = 2 ; RETURN
      END IF

!  there are non-trivial group functions

      DO i = 1, MAX( data%nel, data%ng )
        data%ICALCF( i ) = i
      END DO

!  evaluate the element function values

      CALL ELFUN( data%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,          &
                  data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,          &
                  data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,           &
                  data%lelvar, data%lntvar, data%lstadh, data%lstep,           &
                  data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,          &
                  1, ifstat )
      IF ( ifstat /= 0 ) GO TO 930

!  evaluate the element function values

      CALL ELFUN( data%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,          &
                  data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,          &
                  data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,           &
                  data%lelvar, data%lntvar, data%lstadh, data%lstep,           &
                  data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,          &
                  3, ifstat )
      IF ( ifstat /= 0 ) GO TO 930

!  compute the group argument values ft

      DO ig = 1, data%ng
        ftt = - data%B( ig )

!  include the contribution from the linear element

        DO j = data%ISTADA( ig ), data%ISTADA( ig + 1 ) - 1
          ftt = ftt + data%A( j ) * X( data%ICNA( j ) )
        END DO

!  include the contributions from the nonlinear elements

        DO j = data%ISTADG( ig ), data%ISTADG( ig + 1 ) - 1
          ftt = ftt + data%ESCALE( j ) * data%FUVALS( data%IELING( j ) )
        END DO
        data%FT( ig ) = ftt

!  record the derivatives of trivial groups

        IF ( data%GXEQX( ig ) ) THEN
          data%GVALS( ig, 2 ) = 1.0_wp
          data%GVALS( ig, 3 ) = 0.0_wp
        END IF
      END DO

!  evaluate the group derivative values

      IF ( .NOT. data%altriv ) THEN
        CALL GROUP( data%GVALS, data%ng, data%FT, data%GPVALU, data%ng,        &
                    data%ITYPEG, data%ISTGP, data%ICALCF, data%ltypeg,         &
                    data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,         &
                   .TRUE., igstat )
        IF ( igstat /= 0 ) GO TO 930
      END IF

!  define the real work space needed for ELGRD. Ensure that there is 
!  sufficient space

      IF ( data%lwk2 < data%ng ) THEN
         IF ( data%out > 0 ) WRITE( data%out, 2000 )
         status = 2 ; RETURN
      END IF

!  for unconstrained problems, skip construction of data%WRK( ig ) and skip 
!  specialized construction of gradient and Jacobian.  Call ELGRD instead

!  change the group weightings to include the contributions from the
!  Lagrange multipliers

      IF ( data%numcon > 0 ) THEN
        DO ig = 1, data%ng
          i = data%KNDOFC( ig )
          IF ( i == 0 ) THEN
            data%GSCALE_used( ig ) = data%GSCALE( ig )
          ELSE
            data%GSCALE_used( ig ) = data%GSCALE( ig ) * Y( i )
          END IF
        END DO

!  compute the gradient values. Initialize the gradient and Jacobian (or its 
!  transpose) as zero

        G( : n ) = 0.0_wp
        IF ( jtrans ) THEN
          CJAC( : n, : m ) = 0.0_wp
        ELSE
          CJAC( : m, : n ) = 0.0_wp
        END IF

!  consider the ig-th group

        DO ig = 1, data%ng
          ig1 = ig + 1
          icon = data%KNDOFC( ig )
          istrgv = data%IWORK( data%lstagv + ig )
          iendgv = data%IWORK( data%lstagv + ig1 ) - 1
          nelow = data%ISTADG( ig )
          nelup = data%ISTADG( ig1 ) - 1
          nontrv = .NOT. data%GXEQX( ig )

!  compute the first derivative of the group

          gi = data%GSCALE( ig )
          gii = data%GSCALE_used( ig )
          IF  ( nontrv ) THEN
            gi = gi * data%GVALS( ig, 2 )
            gii = gii * data%GVALS( ig, 2 )
          END IF

!  this is the first gradient evaluation or the group has nonlinear elements

          IF ( data%firstg .OR. nelow <= nelup ) THEN
            data%WRK( data%IWORK( data%lsvgrp + istrgv :                       &
                                  data%lsvgrp + iendgv ) ) = 0.0_wp

!  loop over the group's nonlinear elements

             DO ii = nelow, nelup
               iel = data%IELING( ii )
               k = data%INTVAR( iel )
               l = data%ISTAEV( iel )
               nvarel = data%ISTAEV( iel + 1 ) - l
               scalee = data%ESCALE( ii )
               IF ( data%INTREP( iel ) ) THEN

!  the iel-th element has an internal representation

                 nin = data%INTVAR( iel + 1 ) - k
                 CALL RANGE ( iel, .TRUE., data%FUVALS( k ),                   &
                              data%WRK( n + 1 ), nvarel, nin,                  &
                              data%ITYPEE( iel ), nin, nvarel )
!DIR$ IVDEP
                 DO i = 1, nvarel
                   j = data%IELVAR( l )
                   data%WRK( j ) = data%WRK( j ) + scalee * data%WRK( n + i )
                   l = l + 1
                 END DO
               ELSE

!  the iel-th element has no internal representation

!DIR$ IVDEP
                 DO i = 1, nvarel
                   j = data%IELVAR( l )
                   data%WRK( j ) = data%WRK( j ) + scalee * data%FUVALS( k )
                   k = k + 1 ; l = l + 1
                 END DO
               END IF
             END DO

!  include the contribution from the linear element

!DIR$ IVDEP
            DO k = data%ISTADA( ig ), data%ISTADA( ig1 ) - 1
              j = data%ICNA( k )
              data%WRK( j ) = data%WRK( j ) + data%A( k )
            END DO

!  Allocate a gradient

!DIR$ IVDEP
             DO i = istrgv, iendgv
               ll = data%IWORK( data%lsvgrp + i )

!  the group belongs to the objective function

               IF ( icon == 0 ) THEN
                 G( ll ) = G( ll ) + gi * data%WRK( ll )

!  the group defines a constraint

               ELSE
                 IF ( JTRANS ) THEN
                   CJAC( ll, icon ) = gi * data%WRK( ll )
                 ELSE
                   CJAC( icon, ll ) = gi * data%WRK( ll )
                 END IF
                 IF ( grlagf ) G( ll ) = G( ll ) + gii * data%WRK( ll )
               END IF

!  if the group is non-trivial, also store the nonzero entries of the
!  gradient of the function in GRJAC

               IF ( nontrv ) THEN
                 jj = data%IWORK( data%lstajc + ll )
                 data%FUVALS( data%lgrjac + jj ) = data%WRK( ll )

!  increment the address for the next nonzero in the column of the Jacobian 
!  for variable ll

                 data%IWORK( data%lstajc + ll ) = jj + 1
               END IF
             END DO

!  This is not the first gradient evaluation and there is only a linear element

          ELSE

!  allocate a gradient

!DIR$ IVDEP
            DO k = data%ISTADA( ig ), data%ISTADA( ig1 ) - 1
              ll = data%ICNA( k )

!  the group belongs to the objective function

              IF ( icon == 0 ) THEN
                G( ll ) = G( ll ) + gi * data%A( k )

!  the group defines a constraint

              ELSE
                IF ( JTRANS ) THEN
                  CJAC( ll, icon ) = gi * data%A( k )
                ELSE
                  CJAC( icon, ll ) = gi * data%A( k )
                END IF
                IF ( GRLAGF ) G( ll ) = G( ll ) + gii * data%A( k )
              END IF
            END DO

!  the group is non-trivial; increment the starting addresses for the groups 
!  used by each variable in the (unchanged) linear element to avoid resetting
!  the nonzeros in the Jacobian

            IF ( nontrv ) THEN
!DIR$ IVDEP
              DO i = istrgv, iendgv
                ll = data%IWORK( data%lsvgrp + i )
                data%IWORK( data%lstajc + ll )                               &
                  = data%IWORK( data%lstajc + ll ) + 1
              END DO
            END IF
          END IF
        END DO

!  reset the starting addresses for the lists of groups using each variable to 
!  their values on entry

        DO i = n, 2, - 1
           data%IWORK( data%lstajc + i ) = data%IWORK( data%lstajc + i - 1 )
        END DO
        data%IWORK( data%lstajc + 1 ) = 1
      ELSE

!  compute the gradient value

        CALL ELGRD( n, data%ng, data%firstg, data%ICNA( 1 ), data%licna, &
                      data%ISTADA( 1 ), data%lstada, data%IELING( 1 ), &
                      data%leling, data%ISTADG( 1 ), data%lstadg, &
                      data%ITYPEE( 1 ), data%lintre, &
                      data%ISTAEV( 1 ), data%lstaev, data%IELVAR( 1 ), &
                      data%lelvar, data%INTVAR( 1 ), data%lntvar, &
                      data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, &
                      data%IWORK( data%lstajc + 1 ), data%lnstjc, &
                      data%IWORK( data%lstagv + 1 ), data%lnstgv, &
                      data%A( 1 ), data%la, &
                      data%GVALS( : , 2 ), data%lgvals, &
                      data%FUVALS, data%lnguvl, &
                      data%FUVALS( data%lggfx + 1 ), &
                      data%GSCALE( 1 ), data%lgscal, &
                      data%ESCALE( 1 ), data%lescal, &
                      data%FUVALS( data%lgrjac + 1 ), &
                      data%lngrjc, data%WRK( 1 ), data%WRK( n + 1 ), &
                      data%maxsel, &
                      data%GXEQX( 1 ), data%lgxeqx, &
                      data%INTREP( 1 ), data%lintre, RANGE )

!  store the gradient value

        DO i = 1, n
          G( i ) = data%FUVALS( data%lggfx + i )
        END DO
      END IF
      data%firstg = .FALSE.

!  define the real work space needed for ASMBL. Ensure that there is 
!  sufficient space

      IF ( data%numcon > 0 ) THEN
        lwkh = data%lwk2 - n - 3 * data%maxsel - data%ng
      ELSE
        lwkh = data%lwk2 - n - 3 * data%maxsel
      END IF
      IF ( lwkh <= 0 ) THEN
        IF ( data%out > 0 ) WRITE( data%out, 2000 )
        status = 2 ; RETURN
      END IF

!  define the integer work space needed for ASMBL. Ensure that there is 
!  sufficient space

      liwkh = data%liwk2 - n
      lh = MIN( lwkh, ( liwkh - 3 * n ) / 4 )
      linxtr = lh + n
      IF ( lh <= 0 ) THEN
        IF ( data%out > 0 ) WRITE( data%out, 2010 )
        status = 2 ; RETURN
      END IF

!  set starting addresses for partitions of the integer workspace

      lih = lh
      lirnh = 0
      ljcnh = lirnh + lih
      lnxtrw = ljcnh + lih
      DO i = 1, n
        data%IVAR( i ) = i
      END DO

!  assemble the Hessian

      IF ( data%numcon > 0 ) THEN
      CALL ASMBL( n, data%ng, data%maxsel, n, lh, lih, nnzh, &
                   n, data%IVAR( 1), data%ISTADH( 1 ), data%lstadh, &
                   data%ICNA( 1 ), data%licna, &
                   data%ISTADA( 1 ), data%lstada, &
                   data%INTVAR( 1 ), data%lntvar, &
                   data%IELVAR( 1 ), data%lelvar, &
                   data%IELING( 1 ), data%leling, &
                   data%ISTADG( 1 ), data%lstadg, &
                   data%ISTAEV( 1 ), data%lstaev, &
                   data%IWORK( data%lstagv + 1 ), data%lnstgv, &
                   data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, &
                   data%IWORK( data%lsend + lirnh + 1 ), &
                   data%IWORK( data%lsend + ljcnh + 1 ), &
                   data%IWORK( data%lsend + lnxtrw + 1 ), linxtr, &
                   data%IWORK( data%lsend + liwkh + 1 ), n, &
                   data%A( 1 ), data%la, &
                   data%FUVALS, data%lnguvl, data%FUVALS, data%lnhuvl, &
                   data%GVALS( : , 2 ), data%GVALS( : , 3 ), &
                   data%WRK( 1 ), data%ESCALE( 1 ), data%lescal, &
                   data%WRK( data%ng + 1 ), data%WRK( lwkh + data%ng + 1 ), &
                   data%lwk2 - lwkh, data%GXEQX( 1 ), &
                   data%lgxeqx, data%INTREP( 1 ), data%lintre, &
                   data%ITYPEE( 1 ), data%lintre, &
                   RANGE, 1, data%out, .FALSE., i, inform, &
                   .FALSE., .FALSE. )
      ELSE
      CALL ASMBL( n, data%ng, data%maxsel, n, lh, lih, nnzh, &
                   n, data%IVAR( 1), data%ISTADH( 1 ), data%lstadh, &
                   data%ICNA( 1 ), data%licna, &
                   data%ISTADA( 1 ), data%lstada, &
                   data%INTVAR( 1 ), data%lntvar, &
                   data%IELVAR( 1 ), data%lelvar, &
                   data%IELING( 1 ), data%leling, &
                   data%ISTADG( 1 ), data%lstadg, &
                   data%ISTAEV( 1 ), data%lstaev, &
                   data%IWORK( data%lstagv + 1 ), data%lnstgv, &
                   data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, &
                   data%IWORK( data%lsend + lirnh + 1 ), &
                   data%IWORK( data%lsend + ljcnh + 1 ), &
                   data%IWORK( data%lsend + lnxtrw + 1 ), linxtr, &
                   data%IWORK( data%lsend + liwkh + 1 ), n, &
                   data%A( 1 ), data%la, &
                   data%FUVALS, data%lnguvl, data%FUVALS, data%lnhuvl, &
                   data%GVALS( : , 2 ), data%GVALS( : , 3 ), &
                   data%GSCALE( 1 ), data%ESCALE( 1 ), data%lescal, &
                   data%WRK( 1 ), data%WRK( lwkh + 1 ), &
                   data%lwk2 - lwkh, data%GXEQX( 1 ), &
                   data%lgxeqx, data%INTREP( 1 ), data%lintre, &
                   data%ITYPEE( 1 ), data%lintre, &
                   RANGE, 1, data%out, .FALSE., i, inform, &
                   .FALSE., .FALSE. )
      END IF

!  check that there is sufficient integer workspace

      IF ( inform > 0 ) THEN
        IF ( data%out > 0 ) WRITE( iout, 2010 )
        status = 2 ; RETURN
      END IF

!  initialize the dense Hessian matrix

      H( : n, : n ) = 0.0_wp

!  transfer the matrix from co-ordinate to dense storage and symmetrize the 
!  martix

      IF ( data%numcon > 0 ) THEN
        DO k = 1, nnzh
          i = data%IWORK( data%lsend + lirnh + k )
          j = data%IWORK( data%lsend + ljcnh + k )
          H( i, j ) = data%WRK( data%ng + k )
          H( j, i ) = data%WRK( data%ng + k )
        END DO
      ELSE
        DO k = 1, nnzh
          i = data%IWORK( data%lsend + lirnh + k )
          j = data%IWORK( data%lsend + ljcnh + k )
          H( i, j ) = data%WRK( k )
          H( j, i ) = data%WRK( k )
        END DO
      END IF

!  update the counters for the report tool

      data%nc2og = data%nc2og + 1
      data%nc2cg = data%nc2cg + data%pnc
      data%nc2oh = data%nc2oh + 1
      data%nc2ch = data%nc2ch + data%pnc
      status = 0
      RETURN

!  unsuccessful returns

  930 CONTINUE
      IF ( data%out > 0 ) WRITE( data%out,                                     &
        "( ' ** SUBROUTINE CGRDH: error flag raised during SIF evaluation' )" )
      status = 3
      RETURN

!  non-executable statements

 2000 FORMAT( ' ** SUBROUTINE CGRDH: Increase the size of WK' )
 2010 FORMAT( ' ** SUBROUTINE CGRDH: Increase the size of IWK' )
 2020 FORMAT( ' ** SUBROUTINE CGRDH: Increase the leading dimension of CJAC' )
 2030 FORMAT( ' ** SUBROUTINE CGRDH: Increase the second dimension of CJAC' )
 2040 FORMAT( ' ** SUBROUTINE CGRDH: Increase the leading dimension of H' )

!  end of subroutine CGRDH

      END SUBROUTINE CGRDH
