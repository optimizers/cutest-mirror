! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE UBANDH( data, status, n, goth, X, nsemib, BANDH, lbandh )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n, nsemib, lbandh
      LOGICAL :: goth
      INTEGER, INTENT( OUT ) :: status
      REAL ( KIND = wp ) :: X( n ), BANDH( 0 : lbandh, n )

!  Compute the portion of the Hessian matrix of a group partially
!  separable function which lies within a band of given semi-bandwidth.
!  The diagonal and subdiagonal entries in column i are stored in
!  locations BAND( j, i ), where j = 0 for the diagonal and j > 0 for
!  the j-th subdiagonal entry, j = 1, ..., MIN( nsemib, n - i ) and
!  nsemib is the requested semi-bandwidth.

!  Based on the minimization subroutine data%laNCELOT/SBMIN
!  by Conn, Gould and Toint.

!  Nick Gould, for CGT productions.
!  August 1993.

!  Local variables

      INTEGER :: i, ig, j, lh2, nsemiw, inform, nnzh
      INTEGER :: ifstat, igstat
      INTEGER :: IRNH(   1 ), ICNH(   1 )
      REAL ( KIND = wp ) :: zero, one, ftt
      PARAMETER ( zero = 0.0_wp, one = 1.0_wp )
      EXTERNAL :: RANGE 

!  Check that there is room for the band matrix.

      nsemiw = MAX( 0, MIN( nsemib, n - 1 ) )
      IF ( lbandh < nsemiw ) THEN
         IF ( data%out > 0 ) WRITE( data%out, 2020 )
         status = 2 ; RETURN
      END IF

!  there are non-trivial group functions.

      DO i = 1, MAX( data%nel, data%ng )
        data%ICALCF( i ) = i
      END DO

!  evaluate the element function values.

      CALL ELFUN( data%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,          &
                  data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,          &
                  data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,           &
                  data%lelvar, data%lntvar, data%lstadh, data%lstep,           &
                  data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,          &
                  1, ifstat )

!  evaluate the element function values.

      CALL ELFUN( data%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,          &
                  data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,          &
                  data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,           &
                  data%lelvar, data%lntvar, data%lstadh, data%lstep,           &
                  data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,          &
                  3, ifstat )

!  compute the group argument values ft.

      DO ig = 1, data%ng
        ftt = - data%B( ig )

!  include the contribution from the linear element.

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
          data%GVALS( ig, 2 ) = one
          data%GVALS( ig, 3 ) = zero
        END IF
      END DO

!  evaluate the group derivative values.

        IF ( .NOT. data%altriv )                                               &
          CALL GROUP( data%GVALS, data%ng, data%FT, data%GPVALU, data%ng,      &
                      data%ITYPEG, data%ISTGP, data%ICALCF, data%ltypeg,       &
                      data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,       &
                      .TRUE., igstat )

!  Compute the gradient value.

      CALL ELGRD( n, data%ng, data%firstg, data%ICNA( 1 ), data%licna, &
                   data%ISTADA( 1 ), data%lstada, data%IELING( 1 ), &
                   data%leling, data%ISTADG( 1 ), data%lstadg, &
                   data%ITYPEE( 1 ), data%lintre, &
                   data%ISTAEV( 1 ), data%lstaev, data%IELVAR( 1 ), &
                   data%lelvar, data%INTVAR( 1 ), data%lntvar, data%IWORK( data%lsvgrp + 1 ), &
                   data%lnvgrp, data%IWORK( data%lstajc + 1 ), data%lnstjc, &
                   data%IWORK( data%lstagv + 1 ), data%lnstgv, data%A( 1 ), data%la, &
                   data%GVALS( : , 2 ), data%lgvals, &
                   data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ), &
                   data%GSCALE( 1 ), data%lgscal, &
                   data%ESCALE( 1 ), data%lescal, data%FUVALS( data%lgrjac + 1 ), &
                   data%lngrjc, data%WRK( 1 ), data%WRK( n + 1 ), data%maxsel, &
                   data%GXEQX( 1 ), data%lgxeqx, &
                   data%INTREP( 1 ), data%lintre, RANGE )
      data%firstg = .FALSE.

!  Define the real work space needed for ASMBL.
!  Ensure that there is sufficient space.

      lh2 = n * ( nsemiw + 1 )
      IF ( data%lwk2 < lh2 + n + 3 * data%maxsel ) THEN
         IF ( data%out > 0 ) WRITE( data%out, 2000 )
         status = 2 ; RETURN
      END IF

!  Set starting addresses for partitions of the integer workspace.

      DO i = 1, n
        data%IVAR( i ) = i
      END DO

!  Assemble the Hessian.

      CALL ASMBL( n, data%ng, data%maxsel, nsemiw, lh2, 1, nnzh, &
                   n, data%IVAR( 1), data%ISTADH( 1 ), data%lstadh, &
                   data%ICNA( 1 ), data%licna, &
                   data%ISTADA( 1 ), data%lstada, data%INTVAR( 1 ), data%lntvar, &
                   data%IELVAR( 1 ), data%lelvar, data%IELING( 1 ), data%leling, &
                   data%ISTADG( 1 ), data%lstadg, data%ISTAEV( 1 ), data%lstaev, &
                   data%IWORK( data%lstagv + 1 ), data%lnstgv, data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, &
                   irnh, ICNH, data%IWORK( data%lsend + 1 ), 1, &
                   data%IWORK( data%liwk2 - n + 1 ), n, &
                   data%A( 1 ), data%la, data%FUVALS, data%lnguvl, data%FUVALS, data%lnhuvl, &
                   data%GVALS( : , 2 ), data%GVALS( : , 3 ), &
                   data%GSCALE( 1 ), data%ESCALE( 1 ), data%lescal, &
                   data%WRK( 1 ), data%WRK( lh2 + 1 ), &
                   data%lwk2 - lh2, data%GXEQX( 1 ), &
                   data%lgxeqx, data%INTREP( 1 ), data%lintre, &
                   data%ITYPEE( 1 ), data%lintre, &
                   RANGE, 1, data%out, .TRUE., i, inform, .FALSE., &
                   .FALSE. )

!  Check that there is sufficient integer workspace.

      IF ( inform > 0 ) THEN
        IF ( data%out > 0 ) WRITE( data%out, 2010 )
        status = 2 ; RETURN
      END IF

!  Place the diagonal entries in their correct positions.

      BANDH( 0, : n ) = data%WRK( : n )

!  Now, place the subdiagonal entries in their correct positions.

      lh2 = n
      DO j = 1, n
        DO i = 1, MIN( nsemiw, n - j )
          BANDH( i, j ) = data%WRK( lh2 + i )
        END DO
        lh2 = lh2 + nsemiw
      END DO
      status = 0
      RETURN

!  Non-executable statements

 2000 FORMAT( ' ** SUBROUTINE UBANDH: Increase the size of WK ' )
 2010 FORMAT( ' ** SUBROUTINE UBANDH: Increase the size of IWK ' )
 2020 FORMAT( ' ** SUBROUTINE UBANDH: array dimension lbandh should be', &
              ' larger than nsemib' )

!  end of UBANDH

      END

