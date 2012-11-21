! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE UGRDH( data, status, n, X, G, lh1, H )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n, lh1
      INTEGER, INTENT( OUT ) :: status
      REAL ( KIND = wp ) :: X( n ), G( n ), H( lh1, n )

!  Compute the gradient and Hessian matrix of a group partially 
!  separable function. The Hessian is stored as a dense symmetric matrix.

!  Based on the minimization subroutine data%laNCELOT/SBMIN
!  by Conn, Gould and Toint.

!  Nick Gould, for CGT productions.
!  December 1990.

!  Local variables

      INTEGER :: i, j, k, lih, nnzh, ifstat, igstat
      INTEGER :: ig, lh, lwkh, liwkh, lirnh, ljcnh, lnxtrw, linxtr, inform
      REAL ( KIND = wp ) :: ftt
      EXTERNAL :: RANGE 

!  Check input parameters.

      IF ( lh1 < n ) THEN
        WRITE( data%out, 2020 )
        status = 2 ; RETURN
      END IF

!  there are non-trivial group functions.

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

!  evaluate the element function values

      CALL ELFUN( data%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,          &
                  data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,          &
                  data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,           &
                  data%lelvar, data%lntvar, data%lstadh, data%lstep,           &
                  data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,          &
                  3, ifstat )

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

      IF ( .NOT. data%altriv )                                                 &
        CALL GROUP( data%GVALS, data%ng, data%FT, data%GPVALU, data%ng,        &
                    data%ITYPEG, data%ISTGP, data%ICALCF, data%ltypeg,         &
                    data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,         &
                    .TRUE., igstat )

!  compute the gradient value

      CALL ELGRD( n, data%ng, data%firstg, data%ICNA( 1 ), data%licna, &
                   data%ISTADA( 1 ), data%lstada, data%IELING( 1 ), &
                   data%leling, data%ISTADG( 1 ), data%lstadg, &
                   data%ITYPEE( 1 ), data%lintre, &
                   data%ISTAEV( 1 ), data%lstaev, data%IELVAR( 1 ), &
                   data%lelvar, data%INTVAR( 1 ), data%lntvar, &
                   data%IWORK( data%lsvgrp + 1 ), &
                   data%lnvgrp, data%IWORK( data%lstajc + 1 ), data%lnstjc, &
                   data%IWORK( data%lstagv + 1 ), data%lnstgv, &
                   data%A( 1 ), data%la, &
                   data%GVALS( : , 2 ), data%lgvals, &
                   data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ), &
                   data%GSCALE( 1 ), data%lgscal, &
                   data%ESCALE( 1 ), data%lescal, &
                   data%FUVALS( data%lgrjac + 1 ), &
                   data%lngrjc, data%WRK( 1 ), data%WRK( n + 1 ), data%maxsel, &
                   data%GXEQX( 1 ), data%lgxeqx, &
                   data%INTREP( 1 ), data%lintre, RANGE )
      data%firstg = .FALSE.

!  define the real work space needed for ASMBL. Ensure that there is 
!  sufficient space

      lwkh = data%lwk2 - n - 3 * data%maxsel
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
                   RANGE, 1, data%out, .FALSE., i, inform, .FALSE., &
                   .FALSE. )

!  check that there is sufficient integer workspace

      IF ( inform > 0 ) THEN
         IF ( data%out > 0 ) WRITE( data%out, 2010 )
         status = 2 ; RETURN
      END IF

!  initialize the dense matrix

      H( : n, : n ) = 0.0_wp

!  transfer the matrix from co-ordinate to dense storage and symmetrize it

      DO k = 1, nnzh
        i = data%IWORK( data%lsend + lirnh + k )
        j = data%IWORK( data%lsend + ljcnh + k )
        H( i, j ) = data%WRK( k )
        H( j, i ) = data%WRK( k )
      END DO

!  store the gradient value

      DO i = 1, n
        G( i ) = data%FUVALS( data%lggfx + i )
      END DO

!  update the counters for the report tool

      data%nc2og = data%nc2og + 1
      data%nc2oh = data%nc2oh + 1
      status = 0
      RETURN

!  non-executable statements

 2000 FORMAT( ' ** SUBROUTINE UGRDH: Increase the size of WK ' )
 2010 FORMAT( ' ** SUBROUTINE UGRDH: Increase the size of IWK ' )
 2020 FORMAT( ' ** SUBROUTINE UGRDH: Increase the leading dimension of H ' )

!  end of UGRDH.

      END
