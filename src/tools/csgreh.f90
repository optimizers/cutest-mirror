! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CSGREH( data, status, n, m, X, GRLAGF, Y, nnzj, lcjac,        &
                         CJAC, INDVAR, INDFUN, ne, IRNHI, lirnhi, le,          &
                         IPRNHI, HI, lhi, IPRHI, BYROWS )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n, m, nnzj, ne, lcjac, le, lirnhi, lhi 
      INTEGER, INTENT( OUT ) :: status
      LOGICAL :: GRLAGF, BYROWS
      INTEGER :: INDVAR( LCJAC), INDFUN( lcjac )
      INTEGER :: IRNHI ( lirnhi )
      INTEGER :: IPRNHI( le ), IPRHI ( le )
      REAL ( KIND = wp ) :: X( n ), Y( m ), HI( lhi ), CJAC( lcjac )

!  Compute the Jacobian matrix for an optimization problem
!  initially written in Standard Input Format (SIF).
!  Also compute the Hessian matrix of the Lagrangian function of
!  the problem.

!  The Jacobian is represented in "co-ordinate" format.
!  The Hessian is represented in "finite element format", i.e., 

!           ne
!      H = sum H_i, 
!          i=1

!  where each element H_i involves a small subset of the rows of H.
!  H is stored as a list of the row indices involved in each element
!  and the upper triangle of H_i (stored by rows or columns). 
!  Specifically,

!  ne (integer) number of elements
!  IRNHI (integer array) a list of the row indices involved which each
!          element. Those for element i directly proceed those for 
!          element i + 1, i = 1, ..., NE-1
!  IPRNHI (integer array) pointers to the position in IRNHI of the first 
!          row index in each element. IPRNHI(NE + 1) points to the first 
!          empty location in IRPNHI
!  HI (real array) a list of the nonzeros in the upper triangle of
!          H_i, stored by rows, or by columns, for each element. Those 
!          for element i directly proceed those for element, i + 1, 
!          i = 1, ..., NE-1
!  IPRHI (integer array) pointers to the position in HI of the first 
!          nonzero in each element. IPRHI(NE + 1) points to the first 
!          empty location in HI
!  BYROWS (logical) must be set .TRUE. if the upper triangle of each H_i 
!          is to be stored by rows, and .FALSE. if it is to be stored
!          by columns.

!  Based on the minimization subroutine data%laNCELOT/SBMIN
!  by Conn, Gould and Toint.

!  Nick Gould, for CGT productions,
!  November 1994.

!  Local variables

      INTEGER :: i, j, iel, k, ig, ii, ig1, l, jj, ll, liwkh, icon, lgtemp
      INTEGER :: nin, nvarel, nelow, nelup, istrgv, inform, iendgv
      INTEGER :: ifstat, igstat
      LOGICAL :: nontrv
      REAL ( KIND = wp ) :: ftt, gi, scalee, gii
      EXTERNAL :: RANGE

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

!  evaluate the element function gradients and Hessians

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

!  include the contributions from the nonlinear elements.

        DO j = data%ISTADG( ig ), data%ISTADG( ig + 1 ) - 1
          ftt = ftt + data%ESCALE( j ) * data%FUVALS( data%IELING( j ) )
        END DO
        data%FT( ig ) = ftt

!  Record the derivatives of trivial groups.

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

!  compute the gradient values. Initialize the gradient of the objective 
!  function as zero

         nnzj = 0
         lgtemp = n + data%maxsel
         DO j = 1, n
           data%WRK( lgtemp + j ) = 0.0_wp
         END DO

!  consider the IG-th group

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
           IF ( nontrv ) THEN
             gi = gi  * data%GVALS( ig, 2 )
             gii = gii * data%GVALS( ig, 2 )
           END IF
           data%WRK( data%IWORK( data%lsvgrp + istrgv :                        &
                                 data%lsvgrp + iendgv ) ) = 0.0_wp

!  This is the first gradient evaluation or the group has nonlinear
!  elements.

           IF ( data%firstg .OR. nelow <= nelup ) THEN

!  loop over the group's nonlinear elements

             DO ii = nelow, nelup
               iel = data%IELING( ii )
               k = data%INTVAR( iel )
               l = data%ISTAEV( iel )
               nvarel = data%ISTAEV( iel + 1 ) - l
               scalee = data%ESCALE( ii )
               IF ( data%INTREP( iel ) ) THEN

!  the IEL-th element has an internal representation

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

!  allocate a gradient

!DIR$ IVDEP
             DO i = istrgv, iendgv
               ll = data%IWORK( data%lsvgrp + i )

!  the group belongs to the objective function

               IF ( icon == 0 ) THEN
                 data%WRK( lgtemp + ll )                                       &
                   = data%WRK( lgtemp + ll ) + gi * data%WRK( ll )

!  the group defines a constraint

               ELSE
                 nnzj = nnzj + 1
                 IF ( nnzj <= lcjac ) THEN
                   CJAC ( nnzj ) = gi * data%WRK( ll )
                   INDFUN( nnzj ) = icon
                   INDVAR( nnzj ) = ll
                 END IF
                 IF ( GRLAGF ) data%WRK( lgtemp + ll )                         &
                   = data%WRK( lgtemp + ll ) + gii * data%WRK( ll )
               END IF

!  if the group is non-trivial, also store the nonzero entries of the
!  gradient of the function in GRJAC.

               IF ( nontrv ) THEN
                 jj = data%IWORK( data%lstajc + ll )
                 data%FUVALS( data%lgrjac + jj ) = data%WRK( ll )

!  increment the address for the next nonzero in the column of
!  the jacobian for variable ll

                 data%IWORK( data%lstajc + ll ) = jj + 1
               END IF
             END DO

!  this is not the first gradient evaluation and there is only a linear element

           ELSE

!  include the contribution from the linear element

!DIR$ IVDEP
             DO k = data%ISTADA( ig ),data%ISTADA( ig1 ) - 1
               j = data%ICNA( k )
               data%WRK( j ) = data%WRK( j ) + data%A( k )
             END DO

!  Allocate a gradient.

!DIR$ IVDEP
             DO i = istrgv, iendgv
               ll = data%IWORK( data%lsvgrp + i )

!  the group belongs to the objective function

               IF ( icon == 0 ) THEN
                 data%WRK( lgtemp + ll )                                       &
                   = data%WRK( lgtemp + ll ) + gi * data%WRK( ll )

!  the group defines a constraint

               ELSE
                 nnzj = nnzj + 1
                 IF ( nnzj <= lcjac ) THEN
                   CJAC ( nnzj ) = gi * data%WRK( ll )
                   INDFUN( nnzj ) = icon
                   INDVAR( nnzj ) = ll
                 END IF
                 IF ( GRLAGF ) data%WRK( lgtemp + ll )                        &
                   = data%WRK( lgtemp + ll ) + gii * data%WRK( ll )
               END IF

!  increment the address for the next nonzero in the column of the Jacobian 
!  for variable ll

               IF ( nontrv ) THEN
                 jj = data%IWORK( data%lstajc + ll )
                 data%IWORK( data%lstajc + ll ) = jj + 1
               END IF
             END DO
           END IF
         END DO

!  reset the starting addresses for the lists of groups using each variable to
!  their values on entry

         DO i = n, 2, - 1
           data%IWORK( data%lstajc + i ) = data%IWORK( data%lstajc + i - 1 )
         END DO
         data%IWORK( data%lstajc + 1 ) = 1

!  transfer the gradient of the objective function to the sparse storage scheme

        DO i = 1, n
          nnzj = nnzj + 1
          IF ( nnzj <= lcjac ) THEN
            CJAC ( nnzj ) = data%WRK( lgtemp + i )
            INDFUN( nnzj ) = 0
            INDVAR( nnzj ) = i
          END IF
        END DO
      ELSE

!  compute the gradient value

        CALL CUTEST_form_gradients( n, data%ng, data%nel, data%ntotel,         &
               data%nvrels, data%nnza, data%nvargp, data%firstg, data%ICNA,    &
               data%ISTADA, data%IELING, data%ISTADG, data%ISTAEV,             &
               data%IELVAR, data%INTVAR, data%A, data%GVALS( : , 2 ),          &
               data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ),        &
               data%GSCALE, data%ESCALE, data%FUVALS( data%lgrjac + 1 ),       &
               data%GXEQX, data%INTREP, data%ISVGRP, data%ISTAGV, data%ITYPEE, &
               data%ISTAJC, data%W_ws, data%W_el, RANGE )

!        CALL ELGRD( n, data%ng, data%firstg, data%ICNA( 1 ), data%licna, &
!                      data%ISTADA( 1 ), data%lstada, data%IELING( 1 ), &
!                      data%leling, data%ISTADG( 1 ), data%lstadg, &
!                      data%ITYPEE( 1 ), data%lintre, &
!                      data%ISTAEV( 1 ), data%lstaev, data%IELVAR( 1 ), &
!                      data%lelvar, data%INTVAR( 1 ), data%lntvar, &
!                      data%IWORK( data%lsvgrp + 1 ), &
!                      data%lnvgrp, data%IWORK( data%lstajc + 1 ), data%lnstjc,&
!                      data%IWORK( data%lstagv + 1 ), data%lnstgv, &
!                      data%A( 1 ), data%la, &
!                      data%GVALS( : , 2 ), data%lgvals, &
!                      data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ),&
!                      data%GSCALE( 1 ), data%lgscal, &
!                      data%ESCALE( 1 ), data%lescal, &
!                      data%FUVALS( data%lgrjac + 1 ), &
!                      data%lngrjc, data%WRK( 1 ), data%WRK( n + 1 ), &
!                      data%maxsel, &
!                      data%GXEQX( 1 ), data%lgxeqx, &
!                      data%INTREP( 1 ), data%lintre, RANGE )

!  transfer the gradient of the objective function to the sparse storage scheme

         nnzj = 0
         DO i = 1, n
           nnzj = nnzj + 1
           IF ( nnzj <= lcjac ) THEN
             CJAC ( nnzj ) = data%FUVALS( data%lggfx + i )
             INDFUN( nnzj ) = 0
             INDVAR( nnzj ) = i
           END IF
         END DO
      END IF
      data%firstg = .FALSE.

!  verify that the Jacobian can fit in the alloted space

      IF ( nnzj > lcjac ) THEN
         IF ( data%out > 0 ) WRITE( data%out, 2030 ) nnzj - lcjac 
         status = 2 ; RETURN
      END IF

!  define the real work space needed for ASMBE. Ensure that there is 
!  sufficient space

      IF ( data%numcon > 0 ) THEN
         IF ( data%lwk2 < n + 3 * data%maxsel + data%ng ) THEN
            IF ( data%out > 0 ) WRITE( data%out, 2000 )
            status = 2 ; RETURN
         END IF
      ELSE
         IF ( data%lwk2 < n + 3 * data%maxsel ) THEN
            IF ( data%out > 0 ) WRITE( data%out, 2000 )
            status = 2 ; RETURN
         END IF
      END IF

!  define the integer work space needed for ASMBE.  Ensure that there is 
!  sufficient space

      liwkh = data%liwk2 - n

!  assemble the Hessian

      IF ( data%numcon > 0 ) THEN
      CALL ASMBE( n, data%ng, data%maxsel,  &
                      data%ISTADH( 1 ), data%lstadh, &
                      data%ICNA( 1 ), data%licna, &
                      data%ISTADA( 1 ), data%lstada, &
                      data%INTVAR( 1 ), data%lntvar, &
                      data%IELVAR( 1 ), data%lelvar, &
                      data%IELING( 1 ), data%leling, &
                      data%ISTADG( 1 ), data%lstadg, &
                      data%ISTAEV( 1 ), data%lstaev, &
                      data%IWORK( data%lstagv + 1 ), data%lnstgv, &
                      data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, &
                      data%IWORK( liwkh + 1 ), data%liwk2 - liwkh, &
                      data%A( 1 ), data%la, data%FUVALS, &
                      data%lnguvl, data%FUVALS, data%lnhuvl, &
                      data%GVALS( : , 2 ), data%GVALS( : , 3 ), &
                      data%WRK( 1 ), data%ESCALE( 1 ), data%lescal, &
                      data%WRK( data%ng + 1 ), data%lwk2 - data%ng, &
                      data%GXEQX( 1 ), data%lgxeqx, data%INTREP( 1 ), &
                      data%lintre, data%ITYPEE( 1 ), data%lintre, RANGE, ne,  &
                      IRNHI, lirnhi, IPRNHI, HI, lhi, IPRHI, &
                      BYROWS, 1, data%out, inform )
      ELSE
      CALL ASMBE( n, data%ng, data%maxsel,  &
                      data%ISTADH( 1 ), data%lstadh, &
                      data%ICNA( 1 ), data%licna, &
                      data%ISTADA( 1 ), data%lstada, &
                      data%INTVAR( 1 ), data%lntvar, &
                      data%IELVAR( 1 ), data%lelvar, &
                      data%IELING( 1 ), data%leling, &
                      data%ISTADG( 1 ), data%lstadg, &
                      data%ISTAEV( 1 ), data%lstaev, &
                      data%IWORK( data%lstagv + 1 ), data%lnstgv, &
                      data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, &
                      data%IWORK( liwkh + 1 ), data%liwk2 - liwkh, &
                      data%A( 1 ), data%la, data%FUVALS, &
                      data%lnguvl, data%FUVALS, data%lnhuvl, &
                      data%GVALS( : , 2 ), data%GVALS( : , 3 ), &
                      data%GSCALE( 1 ), data%ESCALE( 1 ), data%lescal, &
                      data%WRK( 1 ), data%lwk2 - data%ng, &
                      data%GXEQX( 1 ), data%lgxeqx, data%INTREP( 1 ), &
                      data%lintre, data%ITYPEE( 1 ), data%lintre, RANGE, ne,  &
                      IRNHI, lirnhi, IPRNHI, HI, lhi, IPRHI, &
                      BYROWS, 1, data%out, inform )
      END IF

!  check that there is room for the elements

      IF ( inform > 0 ) THEN
        IF ( data%out > 0 ) WRITE( data%out, 2020 )
        status = 2 ; RETURN
      END IF

!  update the counters for the report tool

      data%nc2cg = data%nc2og + 1
      data%nc2oh = data%nc2oh + 1
      data%nc2cg = data%nc2cg + data%pnc
      data%nc2ch = data%nc2ch + data%pnc
      status = 0
      RETURN

!  unsuccessful returns

  930 CONTINUE
      IF ( data%out > 0 ) WRITE( data%out,                                     &
        "( ' ** SUBROUTINE CSGREH: error flag raised during SIF evaluation' )" )
      status = 3
      RETURN

!  non-executable statements

 2000 FORMAT( ' ** SUBROUTINE CSGREH: Increase the size of WK' )
 2020 FORMAT( ' ** SUBROUTINE CSGREH: Increase the size of',                   &
              ' IPNRHI, IPRHI, IRNHI or HI' )
 2030 FORMAT( /, ' ** SUBROUTINE CSGREH: array length lcjac too small.',       &
              /, ' -- Increase the parameter lcjac by at least ', I0,          &
                 ' and restart.' )

!  end of subroutine CSGREH

      END SUBROUTINE CSGREH
