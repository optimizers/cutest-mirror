! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE UDH ( data, N, X, LH1, H )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: N, LH1
      REAL ( KIND = wp ) :: X( N ), H( LH1, N )

!  Compute the Hessian matrix of a group partially separable function.
!  The Hessian is stored as a dense symmetric matrix.

!  Based on the minimization subroutine data%laNCELOT/SBMIN
!  by Conn, Gould and Toint.

!  Nick Gould, for CGT productions.
!  December 1990.


! ---------------------------------------------------------------------




! ---------------------------------------------------------------------



! ---------------------------------------------------------------------


! ---------------------------------------------------------------------

!  integer variables from the GLOBAL common block.


!  integer variables from the LOCAL common block.


!  integer variables from the PRFCTS common block.


!  Local variables

      INTEGER :: I, J, K, LIH, NNZH, LJCNH, IFSTAT
      INTEGER :: IG, LH, LWKH, LIWKH, LIRNH, IGSTAT
      INTEGER :: LNXTRW, LINXTR, INFORM
      REAL ( KIND = wp ) :: ZERO, ONE, FTT
      PARAMETER ( ZERO = 0.0_wp, ONE = 1.0_wp )
      EXTERNAL :: RANGE 

!  Check input parameters.

      IF ( LH1 < N ) THEN
         WRITE( IOUT, 2020 )
         STOP
      END IF

!  there are non-trivial group functions.

      DO 10 I = 1, MAX( data%nelnum, data%ng )
        data%ICALCF( I ) = I
   10 CONTINUE

!  evaluate the element function values.

      CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), data%nelnum, &
                   data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                   data%IELVAR( 1 ), data%INTVAR( 1 ), &
                   data%ISTADH( 1 ), data%ISTEP( 1 ), &
                   data%ICALCF( 1 ),  &
                   data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                   data%lntvar, data%lintre, LFUVAL, data%lvscal, data%lepvlu,  &
                   1, IFSTAT )

!  evaluate the element function values.

      CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), data%nelnum, &
                   data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                   data%IELVAR( 1 ), data%INTVAR( 1 ), &
                   data%ISTADH( 1 ), data%ISTEP( 1 ), &
                   data%ICALCF( 1 ),  &
                   data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                   data%lntvar, data%lintre, LFUVAL, data%lvscal, data%lepvlu,  &
                   3, IFSTAT )

!  compute the group argument values ft.

      DO 90 IG = 1, data%ng
         FTT = - data%B( IG )

!  include the contribution from the linear element.

         DO 30 J = data%ISTADA( IG ), data%ISTADA( IG + 1 ) - 1
            FTT = FTT + data%A( J ) * X( data%ICNA( J ) )
   30    CONTINUE

!  include the contributions from the nonlinear elements.

         DO 60 J = data%ISTADG( IG ), data%ISTADG( IG + 1 ) - 1
            FTT = FTT + data%ESCALE( J ) * data%FUVALS( data%IELING( J ) )
   60    CONTINUE
         data%FT( IG ) = FTT

!  Record the derivatives of trivial groups.

         IF ( data%GXEQX( IG ) ) THEN
            data%GVALS( data%ng + IG ) = ONE
            data%GVALS( 2 * data%ng + IG ) = ZERO
         END IF
   90 CONTINUE

!  evaluate the group derivative values.

      IF ( .NOT. data%altriv ) CALL GROUP ( data%GVALS( 1 ), data%ng, &
            data%FT( 1 ), data%GPVALU( 1 ), data%ng, &
            data%ITYPEG( 1 ), data%ISTGP( 1 ), &
            data%ICALCF( 1 ), &
            data%lcalcg, data%ng1, data%lcalcg, data%lcalcg, data%lgpvlu, &
            .TRUE., IGSTAT )

!  Compute the gradient value.

      CALL DELGRD( N, data%ng, data%firstg, data%ICNA( 1 ), data%licna, &
                   data%ISTADA( 1 ), data%lstada, data%IELING( 1 ), &
                   data%leling, data%ISTADG( 1 ), data%lstadg, &
                   data%ITYPEE( 1 ), data%lintre, &
                   data%ISTAEV( 1 ), data%lstaev, data%IELVAR( 1 ), &
                   data%lelvar, data%INTVAR( 1 ), data%lntvar, data%IWORK( data%lsvgrp + 1 ), &
                   data%lnvgrp, data%IWORK( data%lstajc + 1 ), data%lnstjc, &
                   data%IWORK( data%lstagv + 1 ), data%lnstgv, data%A( 1 ), data%la, &
                   data%GVALS( data%ng + 1 ), data%lgvals, &
                   data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ), &
                   data%GSCALE( 1 ), data%lgscal, &
                   data%ESCALE( 1 ), data%lescal, data%FUVALS( data%lgrjac + 1 ), &
                   data%lngrjc, data%WRK( 1 ), data%WRK( N + 1 ), data%maxsel, &
                   data%GXEQX( 1 ), data%lgxeqx, &
                   data%INTREP( 1 ), data%lintre, RANGE )
      data%firstg = .FALSE.

!  Define the real work space needed for ASMBL.
!  Ensure that there is sufficient space.

      LWKH = data%lwk2 - N - 3 * data%maxsel
      IF ( LWKH <= 0 ) THEN
         IF ( IOUT > 0 ) WRITE( IOUT, 2000 )
         STOP
      END IF

!  Define the integer work space needed for ASMBL.
!  Ensure that there is sufficient space.

      LIWKH = data%liwk2 - N
      LH = MIN( LWKH, ( LIWKH - 3 * N ) / 4 )
      LINXTR = LH + N
      IF ( LH <= 0 ) THEN
         IF ( IOUT > 0 ) WRITE( IOUT, 2010 )
         STOP
      END IF

!  Set starting addresses for partitions of the integer workspace.

      LIH = LH
      LIRNH = 0
      LJCNH = LIRNH + LIH
      LNXTRW = LJCNH + LIH
      DO 100 I = 1, N
         data%IVAR( I ) = I
  100 CONTINUE

!  Assemble the Hessian.

      CALL DASMBL( N, data%ng, data%maxsel, N, LH, LIH, NNZH, &
                   N, data%IVAR( 1), data%ISTADH( 1 ), data%lstadh, &
                   data%ICNA( 1 ), data%licna, &
                   data%ISTADA( 1 ), data%lstada, data%INTVAR( 1 ), data%lntvar, &
                   data%IELVAR( 1 ), data%lelvar, data%IELING( 1 ), data%leling, &
                   data%ISTADG( 1 ), data%lstadg, data%ISTAEV( 1 ), data%lstaev, &
                   data%IWORK( data%lstagv + 1 ), data%lnstgv, data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, &
                   data%IWORK( data%lsend + LIRNH + 1 ), data%IWORK( data%lsend + LJCNH + 1 ), &
                   data%IWORK( data%lsend + LNXTRW + 1 ), LINXTR, &
                   data%IWORK( data%lsend + LIWKH + 1 ), N, &
                   data%A( 1 ), data%la, data%FUVALS, data%lnguvl, data%FUVALS, data%lnhuvl, &
                   data%GVALS( data%ng + 1 ), data%GVALS( 2 * data%ng + 1 ), &
                   data%GSCALE( 1 ), data%ESCALE( 1 ), data%lescal, &
                   data%WRK( 1 ), data%WRK( LWKH + 1 ), &
                   data%lwk2 - LWKH, data%GXEQX( 1 ), &
                   data%lgxeqx, data%INTREP( 1 ), data%lintre, &
                   data%ITYPEE( 1 ), data%lintre, &
                   RANGE, 1, IOUT, .FALSE., I, INFORM, .FALSE., &
                   .FALSE. )

!  Check that there is sufficient integer workspace.

      IF ( INFORM > 0 ) THEN
         IF ( IOUT > 0 ) WRITE( IOUT, 2010 )
         STOP
      END IF

!  Initialize the dense matrix.

      DO 120 J = 1, N
         DO 110 I = 1, N
            H( I, J ) = ZERO
  110    CONTINUE
  120 CONTINUE

!  Transfer the matrix from co-ordinate to dense storage and
!  symmetrize the martix.

      DO 130 K = 1, NNZH
         I = data%IWORK( data%lsend + LIRNH + K )
         J = data%IWORK( data%lsend + LJCNH + K )
         H( I, J ) = data%WRK( K )
         H( J, I ) = data%WRK( K )
  130 CONTINUE

!  Update the counters for the report tool.

      data%nc2oh = data%nc2oh + 1
      RETURN

! Non-executable statements.

 2000 FORMAT( ' ** SUBROUTINE UDH: Increase the size of WK ' )
 2010 FORMAT( ' ** SUBROUTINE UDH: Increase the size of IWK ' )
 2020 FORMAT( ' ** SUBROUTINE UDH: Increase the leading dimension', &
              ' of H ' )

!  end of UDH.

      END
