! (  Last modified on 23 Dec 2000 at 22:01:38  )
      SUBROUTINE DCOPY( n, DX, incx, DY, INCY )

!     COPIES A VECTOR, X, TO A VECTOR, Y.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.

      REAL (  KIND = wp  ) :: DX( * ), DY( * )
      INTEGER :: i, incx, INCY, ix, iy, m, mp1, n

      IF( N<=0 )RETURN
      IF( INCX==1.AND.INCY==1 )GO TO 20

!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!          NOT EQUAL TO 1

      ix = 1
      iy = 1
      IF ( incx < 0 ) ix = ( - n + 1 ) * incx + 1
      IF ( INCY < 0 ) iy = ( - n + 1 ) * INCY + 1
      DO 10 i = 1,N
        DY( iy ) = DX( ix )
        ix = ix + incx
        iy = iy + INCY
   10 CONTINUE
      RETURN

!        CODE FOR BOTH INCREMENTS EQUAL TO 1


!        CLEAN-UP LOOP

   20 m = MOD( n, 7 )
      IF(  m == 0  ) GO TO 40
      DO 30 i = 1, m
        DY( i ) = DX( i )
   30 CONTINUE
      IF(  n < 7  ) RETURN
   40 mp1 = m + 1
      DO 50 i = mp1, n, 7
        DY( i ) = DX( i )
        DY( i + 1 ) = DX( i + 1 )
        DY( i + 2 ) = DX( i + 2 )
        DY( i + 3 ) = DX( i + 3 )
        DY( i + 4 ) = DX( i + 4 )
        DY( i + 5 ) = DX( i + 5 )
        DY( i + 6 ) = DX( i + 6 )
   50 CONTINUE
      RETURN
      END



      REAL (  KIND = wp  ) :: FUNCTION DDOT(  n, DX, incx, DY, INCY  )

!     FORMS THE DOT PRODUCT OF TWO VECTORS.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.

      INTEGER :: i, incx, INCY, ix, iy, m, mp1, n
      REAL (  KIND = wp  ) :: DX( * ), DY( * ), DTEMP, zero
      PARAMETER (  zero = 0.0_wp  )

      DDOT = zero
      DTEMP = zero
      IF( n <= 0 )RETURN
      IF( incx == 1 .AND. INCY == 1 )GO TO 20

!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!          NOT EQUAL TO 1

      ix = 1
      iy = 1
      IF ( incx < 0 )IX = ( - n + 1 ) * incx + 1
      IF ( INCY < 0 )IY = ( - n + 1 ) * INCY + 1
      DO 10 i = 1, n
        DTEMP = DTEMP + DX( ix ) * DY( iy )
        ix = ix + incx
        iy = iy + INCY
   10 CONTINUE
      DDOT = DTEMP
      RETURN

!        CODE FOR BOTH INCREMENTS EQUAL TO 1


!        CLEAN-UP LOOP

   20 m = MOD( n, 5 )
      IF(  m == 0  ) GO TO 40
      DO 30 i = 1, m
        DTEMP = DTEMP + DX( i ) * DY( i )
   30 CONTINUE
      IF( n < 5 ) GO TO 60
   40 mp1 = m + 1
      DO 50 i = mp1, n, 5
        DTEMP = DTEMP + DX( i ) * DY( i ) + DX( i + 1 ) * DY( i + 1 ) +        &
           DX( i + 2 ) * DY( i + 2 ) + DX( i + 3 ) * DY( i + 3 ) +             &
           DX( i + 4 ) * DY( i + 4 )
   50 CONTINUE
   60 CONTINUE
      DDOT = DTEMP
      RETURN
      END



      REAL (  KIND = wp  ) :: FUNCTION DNRM2 ( n, DX, incx )
      INTEGER :: n, incx, next
      REAL (  KIND = wp  ) :: DX(  *  ), CUTLO, CUTHI, HITEST, SUM,            &
                       XMAX, zero, one
      INTRINSIC        ABS, SQRT
      INTEGER :: i, j, nn
      PARAMETER (  zero = 0.0_wp, one = 1.0_wp  )

!     EUCLIDEAN NORM OF THE N-VECTOR STORED IN DX(  ) WITH STORAGE
!     INCREMENT incx .
!     IF    n <= 0 RETURN WITH RESULT = 0.
!     IF n >= 1 THEN incx MUST BE >= 1

!           C.L.LAWSON, 1978 JAN 08

!     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE
!     HOPEFULLY APPLICABLE TO ALL MACHINES.
!         CUTLO = MAXIMUM OF  DSQRT( U/EPS )  OVER ALL KNOWN MACHINES.
!         CUTHI = MINIMUM OF  DSQRT( V )      OVER ALL KNOWN MACHINES.
!     WHERE
!         EPS = SMALLEST NO. SUCH THAT EPS + 1. > 1.
!         U = SMALLEST POSITIVE NO. ( UNDERFLOW LIMIT )
!         V = LARGEST  NO. ( OVERFLOW  LIMIT )

!     BRIEF OUTLINE OF ALGORITHM..

!     PHASE 1    SCANS zero COMPONENTS.
!     MOVE TO PHASE 2 WHEN A COMPONENT is NONZERO AND <= CUTLO
!     MOVE TO PHASE 3 WHEN A COMPONENT is > CUTLO
!     MOVE TO PHASE 4 WHEN A COMPONENT is >= CUTHI/M
!     WHERE m = n FOR X(  ) REAL AND m = 2*N FOR COMPLEX.

!     VALUES FOR CUTLO AND CUTHI..
!     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER
!     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..
!     CUTLO, S.P.   U/EPS = 2**( -102 ) FOR  HONEYWELL.  CLOSE SECONDS ARE
!                   UNIVAC AND DEC AT 2**( -103 )
!                   THUS CUTLO = 2**( -51 ) = 4.44089E-16
!     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.
!                   THUS CUTHI = 2**( 63.5 ) = 1.30438E19
!     CUTLO, D.P.   U/EPS = 2**( -67 ) FOR HONEYWELL AND DEC.
!                   THUS CUTLO = 2**( -33.5 ) = 8.23181D-11
!     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19
      DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /

      IF (  n <= 0 ) THEN
         DNRM2 = zero
      ELSE
         next = 1
         SUM = zero
         nn = n * incx

!  BEGIN MAIN LOOP

         i = 1
   20    CONTINUE
         GO TO ( 30, 50, 70, 110 ), next
   30    CONTINUE
         IF(  ABS(  DX(  i  )  ) > CUTLO  ) GO TO 85
         next = 2
         XMAX = zero

!  PHASE 1.  SUM is zero

   50    CONTINUE
         IF ( DX( i ) == zero  ) GO TO 200
         IF ( ABS( DX( i ) ) > CUTLO  ) GO TO 85

!  PREPARE FOR PHASE 2.

         next = 3
         GO TO 105

!  PREPARE FOR PHASE 4.

  100    CONTINUE
         i = j
         next = 4
         SUM = (  SUM / DX(  i  )  ) / DX(  i  )
  105    CONTINUE
         XMAX = ABS(  DX(  i  )  )
         SUM = SUM + (  DX(  i  ) / XMAX  ) ** 2
         GO TO 200

!  PHASE 2.  SUM is SMALL. SCALE TO AVOID DESTRUCTIVE UNDERFLOW.

   70    CONTINUE
         IF (  ABS(  DX(  i  )  ) > CUTLO  ) THEN

!  PREPARE FOR PHASE 3.

            SUM = (  SUM * XMAX ) * XMAX
            GO TO 85
         END IF

!  COMMON CODE FOR PHASES 2 AND 4.
!  IN PHASE 4 SUM is LARGE.  SCALE TO AVOID OVERFLOW.

  110    CONTINUE
         IF (  ABS(  DX(  i  )  ) > XMAX  ) THEN
            SUM = one + SUM * (  XMAX / DX(  i  )  ) ** 2
            XMAX = ABS(  DX(  i  )  )
         ELSE
            SUM = SUM + (  DX(  i  ) / XMAX  ) ** 2
         END IF
  200    CONTINUE
         i = i + incx
         IF (  i <= nn  ) GO TO 20

!  END OF MAIN LOOP.

!  COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.

         DNRM2 = XMAX * SQRT( SUM )
         GO TO 300

!  FOR REAL OR D.P. SET HITEST = CUTHI/N

   85    CONTINUE
         HITEST = CUTHI / FLOAT(  n  )

!  PHASE 3. SUM is MID-RANGE.  NO SCALING.

         DO 95 j = i, nn, incx
            IF(  ABS(  DX(  j  )  ) >= HITEST  ) GO TO 100
            SUM = SUM + DX(  j  ) ** 2
   95    CONTINUE
         DNRM2 = SQRT(  SUM  )
      END IF
  300 CONTINUE
      RETURN
      END



      SUBROUTINE DAXPY(  n, DA, DX, incx, DY, INCY  )

!     CONSTANT TIMES A VECTOR PLUS A VECTOR.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.

      INTEGER :: i, incx, INCY, ix, iy, m, mp1, n
      REAL (  KIND = wp  ) :: DX( * ), DY( * ), DA, zero

!  SET CONSTANT REAL PARAMETERS.

      PARAMETER (  zero = 0.0_wp  )

      IF( N<=0 )RETURN
      IF ( DA == zero ) RETURN
      IF( INCX==1.AND.INCY==1 )GO TO 20

!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!          NOT EQUAL TO 1

      ix = 1
      iy = 1
      IF( incx < 0 ) ix = ( - n + 1 ) * incx + 1
      IF( INCY < 0 ) iy = ( - n + 1 ) * INCY + 1
      DO 10 i = 1, n
        DY( iy ) = DY( iy ) + DA * DX( ix )
        ix = ix + incx
        iy = iy + INCY
   10 CONTINUE
      RETURN

!        CODE FOR BOTH INCREMENTS EQUAL TO 1


!        CLEAN-UP LOOP

   20 m = MOD( n, 4 )
      IF(  m == 0  ) GO TO 40
      DO 30 i = 1,M
        DY( i ) = DY( i ) + DA * DX( i )
   30 CONTINUE
      IF(  n < 4  ) RETURN
   40 mp1 = m + 1
      DO 50 i = mp1, n, 4
        DY( i ) = DY( i ) + DA * DX( i )
        DY( i + 1 ) = DY( i + 1 ) + DA * DX( i + 1 )
        DY( i + 2 ) = DY( i + 2 ) + DA * DX( i + 2 )
        DY( i + 3 ) = DY( i + 3 ) + DA * DX( i + 3 )
   50 CONTINUE
      RETURN
      END



      SUBROUTINE DROT (  n, DX, incx, DY, INCY, C, S  )

!     APPLIES A PLANE ROTATION.
!     JACK DONGARRA, LINPACK, 3/11/78.

      REAL (  KIND = wp  ) :: DX( * ), DY( * ), DTEMP, C, S
      INTEGER :: i, incx, INCY, ix, iy, n

      IF( n <= 0 )RETURN
      IF( incx == 1 .AND. INCY == 1 )GO TO 20

!       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
!         TO 1

      ix = 1
      iy = 1
      IF( incx < 0 ) ix = ( - n + 1 ) * incx + 1
      IF( INCY < 0 ) iy = ( - n + 1 ) * INCY + 1
      DO 10 i = 1,N
        DTEMP = C * DX( ix ) + S * DY( iy )
        DY( iy ) = C * DY( iy ) - S * DX( ix )
        DX( ix ) = DTEMP
        ix = ix + incx
        iy = iy + INCY
   10 CONTINUE
      RETURN

!       CODE FOR BOTH INCREMENTS EQUAL TO 1

   20 DO 30 i = 1, n
        DTEMP = C * DX( i ) + S * DY( i )
        DY( i ) = C * DY( i ) - S * DX( i )
        DX( i ) = DTEMP
   30 CONTINUE
      RETURN
      END
      SUBROUTINE DROTG( DA, DB, C, S )

!     CONSTRUCT GIVENS PLANE ROTATION.
!     JACK DONGARRA, LINPACK, 3/11/78.

      REAL (  KIND = wp  ) :: DA, DB, C, S, ROE, SCALE, R, Z, zero, one
      INTRINSIC ABS, SQRT, SIGN
      PARAMETER (  zero = 0.0_wp, one = 1.0_wp  )

      ROE = DB
      IF(  ABS( DA ) > ABS( DB )  ) ROE = DA
      SCALE = ABS( DA ) + ABS( DB )
      IF(  SCALE /= zero  ) GO TO 10
         C = one
         S = zero
         R = zero
         GO TO 20
   10 R = SCALE*SQRT( ( DA/SCALE )**2 + ( DB/SCALE )**2 )
      R = SIGN( one,ROE )*R
      C = DA / R
      S = DB / R
   20 Z = one
      IF ( ABS( DA ) > ABS( DB )  ) Z = S
      IF ( ABS( DB ) >= ABS( DA ) .AND. C /= zero  ) Z = one / C
      DA = R
      DB = Z
      RETURN
      END



      SUBROUTINE DSCAL( n, DA, DX, incx )

!     SCALES A VECTOR BY A CONSTANT.
!     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.

      REAL (  KIND = wp  ) :: DA, DX( * )
      INTEGER :: i, incx, m, mp1, n, nincx

      IF( n <= 0 )RETURN
      IF( incx == 1 )GO TO 20

!        CODE FOR INCREMENT NOT EQUAL TO 1

      nincx = n * incx
      DO 10 i = 1, nincx, incx
        DX( i ) = DA * DX( i )
   10 CONTINUE
      RETURN

!        CODE FOR INCREMENT EQUAL TO 1


!        CLEAN-UP LOOP

   20 m = MOD( n, 5 )
      IF(  m == 0  ) GO TO 40
      DO 30 i = 1, m
        DX( i ) = DA * DX( i )
   30 CONTINUE
      IF(  n < 5  ) RETURN
   40 mp1 = m + 1
      DO 50 i = mp1, n, 5
        DX( i ) = DA * DX( i )
        DX( i + 1 ) = DA * DX( i + 1 )
        DX( i + 2 ) = DA * DX( i + 2 )
        DX( i + 3 ) = DA * DX( i + 3 )
        DX( i + 4 ) = DA * DX( i + 4 )
   50 CONTINUE
      RETURN
      END



      SUBROUTINE  DSWAP( n, DX, incx, DY, INCY )

!     INTERCHANGES TWO VECTORS.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.

      REAL (  KIND = wp  ) :: DX( * ), DY( * ), DTEMP
      INTEGER :: i, incx, INCY, ix, iy, m, mp1, n

      IF( N<=0 )RETURN
      IF( INCX==1.AND.INCY==1 )GO TO 20

!       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
!         TO 1

      ix = 1
      iy = 1
      IF( incx < 0 ) ix = ( - n + 1 ) * incx + 1
      IF( INCY < 0 ) iy = ( - n + 1 ) * INCY + 1
      DO 10 i = 1, n
        DTEMP = DX( ix )
        DX( ix ) = DY( iy )
        DY( iy ) = DTEMP
        ix = ix + incx
        iy = iy + INCY
   10 CONTINUE
      RETURN

!       CODE FOR BOTH INCREMENTS EQUAL TO 1


!       CLEAN-UP LOOP

   20 m = MOD( n, 3 )
      IF(  m == 0  ) GO TO 40
      DO 30 i = 1, m
        DTEMP = DX( i )
        DX( i ) = DY( i )
        DY( i ) = DTEMP
   30 CONTINUE
      IF(  n < 3  ) RETURN
   40 mp1 = m + 1
      DO 50 i = mp1, n, 3
        DTEMP = DX( i )
        DX( i ) = DY( i )
        DY( i ) = DTEMP
        DTEMP = DX( i + 1 )
        DX( i + 1 ) = DY( i + 1 )
        DY( i + 1 ) = DTEMP
        DTEMP = DX( i + 2 )
        DX( i + 2 ) = DY( i + 2 )
        DY( i + 2 ) = DTEMP
   50 CONTINUE
      RETURN
      END



      INTEGER FUNCTION IDAMAX( n, DX, incx )

!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array( 1 ) declarations changed to array( * )

      REAL (  KIND = wp  ) :: DX( * ), DMAX
      INTEGER :: i, incx, ix, n

      IDAMAX = 0
      IF(  n < 1 .OR. incx <= 0  ) RETURN
      IDAMAX = 1
      IF( N==1 )RETURN
      IF( INCX==1 )GO TO 20

!        code for increment not equal to 1

      ix = 1
      DMAX = ABS( DX( 1 ) )
      ix = ix + incx
      DO 10 i = 2,N
         IF( ABS( DX( ix ) ) <= DMAX ) GO TO 5
         IDAMAX = i
         DMAX = ABS( DX( ix ) )
    5    ix = ix + incx
   10 CONTINUE
      RETURN

!        code for increment equal to 1

   20 DMAX = ABS( DX( 1 ) )
      DO 30 i = 2,N
         IF( ABS( DX( i ) ) <= DMAX ) GO TO 30
         IDAMAX = i
         DMAX = ABS( DX( i ) )
   30 CONTINUE
      RETURN
      END



      REAL (  KIND = wp  ) :: FUNCTION DASUM(  n, DX, incx  )

!     takes the sum of the absolute values.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array( 1 ) declarations changed to array( * )

      REAL (  KIND = wp  ) :: DX( * ),DTEMP
      INTEGER :: i,INCX,M,MP1,N,NINCX

      DASUM = 0.0D0
      DTEMP = 0.0D0
      IF(  n <= 0 .OR. incx <= 0  )RETURN
      IF( INCX==1 )GO TO 20

!        code for increment not equal to 1

      nincx = N*INCX
      DO 10 i = 1,NINCX,INCX
        DTEMP = DTEMP + ABS( DX( i ) )
   10 CONTINUE
      DASUM = DTEMP
      RETURN

!        code for increment equal to 1


!        clean-up loop

   20 m = MOD( n, 6 )
      IF(  m == 0  ) GO TO 40
      DO 30 i = 1,M
        DTEMP = DTEMP + ABS( DX( i ) )
   30 CONTINUE
      IF(  n < 6  ) GO TO 60
   40 mp1 = m + 1
      DO 50 i = mp1, n, 6
        DTEMP = DTEMP + ABS( DX( i ) ) + ABS( DX( i + 1 ) ) +                  &
         ABS( DX( i + 2 ) ) + ABS( DX( i + 3 ) ) + ABS( DX( i + 4 ) ) +        &
         ABS( DX( i + 5 ) )
   50 CONTINUE
   60 CONTINUE
      DASUM = DTEMP
      RETURN
      END
