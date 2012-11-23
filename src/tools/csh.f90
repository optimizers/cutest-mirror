! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CSH( data, status, n, m, X, Y, nnzh, lh, H, IRNH, ICNH )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n, m, nnzh, lh
      INTEGER, INTENT( OUT ) :: status
      INTEGER :: IRNH( lh ), ICNH( lh )
      REAL ( KIND = wp ) :: X ( n ), Y( m ), H( lh )

!  Compute the Hessian matrix of the Lagrangian function of
!  a problem initially written in Standard Input Format (SIF).

!  The upper triangle of the Hessian is stored in coordinate form,
!  i.e., the entry H(i) has row index IRNH(i) and column index ICNH(i)
!  for i = 1, ...., nnzh.

!  Based on the minimization subroutine data%laNCELOT/SBMIN
!  by Conn, Gould and Toint.

!  Nick Gould, for CGT productions,
!  November 1991.

!  Local variables

      INTEGER :: i, j, ifstat, igstat
      INTEGER :: ig, lih, liwkh, lnxtrw, linxtr, inform
      REAL ( KIND = wp ) :: ftt
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

!  compute the gradient value

        CALL CUTEST_form_gradients( n, data%ng, data%nel, data%ntotel,         &
               data%nvrels, data%nnza, data%nvargp, data%firstg, data%ICNA,    &
               data%ISTADA, data%IELING, data%ISTADG, data%ISTAEV,             &
               data%IELVAR, data%INTVAR, data%A, data%GVALS( : , 2 ),          &
               data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ),        &
               data%GSCALE_used, data%ESCALE, data%FUVALS( data%lgrjac + 1 ),  &
               data%GXEQX, data%INTREP, data%ISVGRP, data%ISTAGV, data%ITYPEE, &
               data%ISTAJC, data%W_ws, data%W_el, RANGE )
      ELSE
        CALL CUTEST_form_gradients( n, data%ng, data%nel, data%ntotel,         &
               data%nvrels, data%nnza, data%nvargp, data%firstg, data%ICNA,    &
               data%ISTADA, data%IELING, data%ISTADG, data%ISTAEV,             &
               data%IELVAR, data%INTVAR, data%A, data%GVALS( : , 2 ),          &
               data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ),        &
               data%GSCALE, data%ESCALE, data%FUVALS( data%lgrjac + 1 ),       &
               data%GXEQX, data%INTREP, data%ISVGRP, data%ISTAGV, data%ITYPEE, &
               data%ISTAJC, data%W_ws, data%W_el, RANGE )
      END IF

!      CALL ELGRD( n, data%ng, data%firstg, data%ICNA( 1 ), data%licna, &
!                  data%ISTADA( 1 ), data%lstada, data%IELING( 1 ), &
!                  data%leling, data%ISTADG( 1 ), data%lstadg, &
!                  data%ITYPEE( 1 ), data%lintre, &
!                  data%ISTAEV( 1 ), data%lstaev, data%IELVAR( 1 ), &
!                  data%lelvar, data%INTVAR( 1 ), data%lntvar, &
!                  data%IWORK( data%lsvgrp + 1 ), &
!                  data%lnvgrp, data%IWORK( data%lstajc + 1 ), data%lnstjc, &
!                  data%IWORK( data%lstagv + 1 ), data%lnstgv, &
!                  data%A( 1 ), data%la, &
!                  data%GVALS( : , 2 ), data%lgvals, &
!                  data%FUVALS, data%lnguvl, &
!                  data%FUVALS( data%lggfx + 1 ), &
!                  data%WRK( 1 ), data%ng, &
!                  data%ESCALE( 1 ), data%lescal, &
!                  data%FUVALS( data%lgrjac + 1 ), &
!                  data%lngrjc, data%WRK( 1 ), data%WRK( n + 1 ), &
!                  data%maxsel, &
!                  data%GXEQX( 1 ), data%lgxeqx, &
!                  data%INTREP( 1 ), data%lintre, RANGE )
!      ELSE
!
!      CALL ELGRD( n, data%ng, data%firstg, data%ICNA( 1 ), data%licna, &
!                  data%ISTADA( 1 ), data%lstada, data%IELING( 1 ), &
!                  data%leling, data%ISTADG( 1 ), data%lstadg, &
!                  data%ITYPEE( 1 ), data%lintre, &
!                  data%ISTAEV( 1 ), data%lstaev, data%IELVAR( 1 ), &
!                  data%lelvar, data%INTVAR( 1 ), data%lntvar, &
!                  data%IWORK( data%lsvgrp + 1 ), &
!                  data%lnvgrp, data%IWORK( data%lstajc + 1 ), data%lnstjc, &
!                  data%IWORK( data%lstagv + 1 ), data%lnstgv, &
!                  data%A( 1 ), data%la, &
!                  data%GVALS( : , 2 ), data%lgvals, &
!                  data%FUVALS, data%lnguvl, &
!                  data%FUVALS( data%lggfx + 1 ), &
!                  data%GSCALE( 1 ), data%lgscal, &
!                  data%ESCALE( 1 ), data%lescal, &
!                  data%FUVALS( data%lgrjac + 1 ), &
!                  data%lngrjc, data%WRK( 1 ), data%WRK( n + 1 ), &
!                  data%maxsel, &
!                  data%GXEQX( 1 ), data%lgxeqx, &
!                  data%INTREP( 1 ), data%lintre, RANGE )

      data%firstg = .FALSE.

!  define the real work space needed for ASMBL. Ensure that there is 
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

!  define the integer work space needed for ASMBL. Ensure that there is 
!  sufficient space

      liwkh = data%liwk2 - n
      linxtr = liwkh / 2
      IF ( linxtr < n ) THEN
         IF ( data%out > 0 ) WRITE( data%out, 2010 )
         status = 2 ; RETURN
      END IF

!  set starting addresses for partitions of the integer workspace

      lih = lh
      lnxtrw = 0
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
                   IRNH, ICNH, data%IWORK( data%lsend + lnxtrw + 1 ), linxtr, &
                   data%IWORK( data%lsend + liwkh + 1 ), n, &
                   data%A( 1 ), data%la, data%FUVALS, data%lnguvl, &
                   data%FUVALS, data%lnhuvl, &
                   data%GVALS( : , 2 ), data%GVALS( : , 3 ), &
                   data%WRK( 1 ), data%ESCALE( 1 ), data%lescal, &
                   H, data%WRK( data%ng + 1 ), data%lwk2 - data%ng, &
                   data%GXEQX( 1 ), data%lgxeqx, data%INTREP( 1 ), &
                   data%lintre, data%ITYPEE( 1 ), data%lintre, &
                   RANGE, 1, data%out, .FALSE., i, inform,  &
                   .FALSE., .TRUE. )
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
                   IRNH, ICNH, data%IWORK( data%lsend + lnxtrw + 1 ), linxtr, &
                   data%IWORK( data%lsend + liwkh + 1 ), n,  &
                   data%A( 1 ), data%la, data%FUVALS, data%lnguvl, &
                   data%FUVALS, data%lnhuvl, &
                   data%GVALS( : , 2 ), data%GVALS( : , 3 ), &
                   data%GSCALE( 1 ), data%ESCALE( 1 ), data%lescal, &
                   H, data%WRK( 1 ), data%lwk2 - data%ng, &
                   data%GXEQX( 1 ), data%lgxeqx, data%INTREP( 1 ), &
                   data%lintre, data%ITYPEE( 1 ), data%lintre, &
                   RANGE, 1, data%out, .FALSE., i, inform, &
                   .FALSE., .TRUE. )
      END IF

!  check that there is sufficient integer workspace

      IF ( inform > 0 ) THEN
         IF ( data%out > 0 ) WRITE( data%out, 2010 )
         status = 2 ; RETURN
      END IF

!  update the counters for the report tool

      data%nc2oh = data%nc2oh + 1
      data%nc2ch = data%nc2ch + data%pnc
      status = 0
      RETURN

!  unsuccessful returns

  930 CONTINUE
      IF ( data%out > 0 ) WRITE( data%out,                                     &
        "( ' ** SUBROUTINE CSH: error flag raised during SIF evaluation' )" )
      status = 3
      RETURN

!  non-executable statements

 2000 FORMAT( ' ** SUBROUTINE CSH: Increase the size of WK ' )
 2010 FORMAT( ' ** SUBROUTINE CSH: Increase the sizes of IWK and LH' )

!  end of subroutine CSH

      END SUBROUTINE CSH
