! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CDH ( data, n, m, X, lv, V, lh1, H )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n, m, lv, lh1
      REAL ( KIND = wp ) :: X( n ), V( lv ), H( lh1, n )

!  Compute the Hessian matrix of the Lagrangian function of
!  a problem initially written in Standard Input Format (SIF).

!  H is a two-dimensional array which gives the value of the
!    Hessian matrix of the Lagrangian function evaluated at
!    X and V. The i,j-th component of the array will contain
!    the derivative with respect to variables X(i) and X(j).

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

      INTEGER :: i, j, k, lih, nnzh, ifstat, igstat
      INTEGER :: ig, lh, lwkh, liwkh, lirnh, ljcnh
      INTEGER :: lnxtrw, linxtr, inform
      REAL ( KIND = wp ) :: zero, one, ftt
      PARAMETER ( zero = 0.0_wp, one = 1.0_wp )
      EXTERNAL :: RANGE 

!  Check input parameters.

      IF ( lh1 < n ) THEN
         WRITE( iout, 2020 )
         STOP
      END IF

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

!  evaluate the element function Hessian values.

      CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), data%nelnum, &
                   data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                   data%IELVAR( 1 ), data%INTVAR( 1 ), &
                   data%ISTADH( 1 ), data%ISTEP( 1 ), &
                   data%ICALCF( 1 ),  &
                   data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                   data%lntvar, data%lintre, lfuval, data%lvscal, data%lepvlu,  &
                   3, ifstat )

!  compute the group argument values ft.

      DO 70 ig = 1, data%ng
         ftt = - data%B( ig )

!  include the contribution from the linear element.

         DO 30 j = data%ISTADA( ig ), data%ISTADA( ig + 1 ) - 1
            ftt = ftt + data%A( j ) * X( data%ICNA( j ) )
   30    CONTINUE

!  include the contributions from the nonlinear elements.

         DO 60 j = data%ISTADG( ig ), data%ISTADG( ig + 1 ) - 1
            ftt = ftt + data%ESCALE( j ) * data%FUVALS( data%IELING( j ) )
   60    CONTINUE
         data%FT( ig ) = ftt

!  Record the derivatives of trivial groups.

         IF ( data%GXEQX( ig ) ) THEN
            data%GVALS( data%ng + ig ) = one
            data%GVALS( 2 * data%ng + ig ) = zero
         END IF
   70 CONTINUE

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
                      data%WRK( 1 ), data%ng, &
                      data%ESCALE( 1 ), data%lescal, data%FUVALS( data%lgrjac + 1 ), &
                      data%lngrjc, data%WRK( 1 ), data%WRK( n + 1 ), data%maxsel, &
                      data%GXEQX( 1 ), data%lgxeqx, &
                      data%INTREP( 1 ), data%lintre, RANGE )
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
      DO 90 i = 1, n
         data%IVAR( i ) = i
   90 CONTINUE

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
                   RANGE, 1, iout, .FALSE., i, inform, .FALSE., &
                   .FALSE. )
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
                   RANGE, 1, iout, .FALSE., i, inform, .FALSE., &
                   .FALSE. )
      END IF

!  Check that there is sufficient integer workspace.

      IF ( inform > 0 ) THEN
         IF ( iout > 0 ) WRITE( iout, 2010 )
         STOP
      END IF

!  Initialize the dense matrix.

      DO 120 j = 1, n
         DO 110 i = 1, n
            H( i, j ) = zero
  110    CONTINUE
  120 CONTINUE

!  Transfer the matrix from co-ordinate to dense storage and
!  symmetrize the martix.

      IF ( data%numcon > 0 ) THEN
         DO 130 k = 1, nnzh
            i = data%IWORK( data%lsend + lirnh + k )
            j = data%IWORK( data%lsend + ljcnh + k )
            H( i, j ) = data%WRK( data%ng + k )
            H( j, i ) = data%WRK( data%ng + k )
  130    CONTINUE
      ELSE
         DO 140 k = 1, nnzh
            i = data%IWORK( data%lsend + lirnh + k )
            j = data%IWORK( data%lsend + ljcnh + k )
            H( i, j ) = data%WRK( k )
            H( j, i ) = data%WRK( k )
  140    CONTINUE
      END IF

!  Update the counters for the report tool.

      data%nc2oh = data%nc2oh + 1
      data%nc2ch = data%nc2ch + data%pnc
      RETURN

! Non-executable statements.

 2000 FORMAT( ' ** SUBROUTINE CDH: Increase the size of WK ' )
 2010 FORMAT( ' ** SUBROUTINE CDH: Increase the size of IWK ' )
 2020 FORMAT( ' ** SUBROUTINE CDH: Increase the leading dimension', &
              ' of H ' )

!  end of CDH.

      END
