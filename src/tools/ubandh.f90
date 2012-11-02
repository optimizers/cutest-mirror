! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE UBANDH( data, N, GOTH, X, NSEMIB, BANDH, LBANDH )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: N, NSEMIB, LBANDH
      LOGICAL :: GOTH
      REAL ( KIND = wp ) :: X( N ), BANDH( 0: LBANDH, N )

!  Compute the portion of the Hessian matrix of a group partially
!  separable function which lies within a band of given semi-bandwidth.
!  The diagonal and subdiagonal entries in column i are stored in
!  locations BAND( j, i ), where j = 0 for the diagonal and j > 0 for
!  the j-th subdiagonal entry, j = 1, ..., MIN( NSEMIB, N - i ) and
!  NSEMIB is the requested semi-bandwidth.

!  Based on the minimization subroutine data%laNCELOT/SBMIN
!  by Conn, Gould and Toint.

!  Nick Gould, for CGT productions.
!  August 1993.


! ---------------------------------------------------------------------




! ---------------------------------------------------------------------



! ---------------------------------------------------------------------


! ---------------------------------------------------------------------

!  integer variables from the GLOBAL common block.


!  integer variables from the LOCAL common block.


!  Local variables

      INTEGER :: I, IG, J, LH2, NSEMIW, INFORM, NNZH
      INTEGER :: IFSTAT, IGSTAT
      INTEGER :: IRNH(   1 ), ICNH(   1 )
      REAL ( KIND = wp ) :: ZERO, ONE, FTT
      PARAMETER ( ZERO = 0.0_wp, ONE = 1.0_wp )
      EXTERNAL :: RANGE 

!  Check that there is room for the band matrix.

      NSEMIW = MAX( 0, MIN( NSEMIB, N - 1 ) )
      IF ( LBANDH < NSEMIW ) THEN
         IF ( IOUT > 0 ) WRITE( IOUT, 2020 )
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

      LH2 = N * ( NSEMIW + 1 )
      IF ( data%lwk2 < LH2 + N + 3 * data%maxsel ) THEN
         IF ( IOUT > 0 ) WRITE( IOUT, 2000 )
         STOP
      END IF

!  Set starting addresses for partitions of the integer workspace.

      DO 100 I = 1, N
         data%IVAR( I ) = I
  100 CONTINUE

!  Assemble the Hessian.

      CALL DASMBL( N, data%ng, data%maxsel, NSEMIW, LH2, 1, NNZH, &
                   N, data%IVAR( 1), data%ISTADH( 1 ), data%lstadh, &
                   data%ICNA( 1 ), data%licna, &
                   data%ISTADA( 1 ), data%lstada, data%INTVAR( 1 ), data%lntvar, &
                   data%IELVAR( 1 ), data%lelvar, data%IELING( 1 ), data%leling, &
                   data%ISTADG( 1 ), data%lstadg, data%ISTAEV( 1 ), data%lstaev, &
                   data%IWORK( data%lstagv + 1 ), data%lnstgv, data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, &
                   IRNH, ICNH, data%IWORK( data%lsend + 1 ), 1, &
                   data%IWORK( data%liwk2 - N + 1 ), N, &
                   data%A( 1 ), data%la, data%FUVALS, data%lnguvl, data%FUVALS, data%lnhuvl, &
                   data%GVALS( data%ng + 1 ), data%GVALS( 2 * data%ng + 1 ), &
                   data%GSCALE( 1 ), data%ESCALE( 1 ), data%lescal, &
                   data%WRK( 1 ), data%WRK( LH2 + 1 ), &
                   data%lwk2 - LH2, data%GXEQX( 1 ), &
                   data%lgxeqx, data%INTREP( 1 ), data%lintre, &
                   data%ITYPEE( 1 ), data%lintre, &
                   RANGE, 1, IOUT, .TRUE., I, INFORM, .FALSE., &
                   .FALSE. )

!  Check that there is sufficient integer workspace.

      IF ( INFORM > 0 ) THEN
         IF ( IOUT > 0 ) WRITE( IOUT, 2010 )
         STOP
      END IF

!  Place the diagonal entries in their correct positions.

      DO 110 J = 1, N
         BANDH( 0, J ) = data%WRK( J )
  110 CONTINUE

!  Now, place the subdiagonal entries in their correct positions.

      LH2 = data%lwkstr + N
      DO 130 J = 1, N
         DO 120 I = 1, MIN( NSEMIW, N - J )
            BANDH( I, J ) = WK( LH2 + I )
  120    CONTINUE
         LH2 = LH2 + NSEMIW
  130 CONTINUE
      RETURN

! Non-executable statements.

 2000 FORMAT( ' ** SUBROUTINE UBANDH: Increase the size of WK ' )
 2010 FORMAT( ' ** SUBROUTINE UBANDH: Increase the size of IWK ' )
 2020 FORMAT( ' ** SUBROUTINE UBANDH: array dimension LBANDH should be', &
              ' larger than NSEMIB' )

!  end of UBANDH.

      END

