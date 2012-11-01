! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)

!     COPIES A VECTOR, X, TO A VECTOR, Y.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.

      REAL ( KIND = wp ) :: DX(*),DY(*)
      INTEGER :: I,INCX,INCY,IX,IY,M,MP1,N

      IF(N<=0)RETURN
      IF(INCX==1.AND.INCY==1)GO TO 20

!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!          NOT EQUAL TO 1

      IX = 1
      IY = 1
      IF(INCX<0)IX = (-N + 1)*INCX + 1
      IF(INCY<0)IY = (-N + 1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN

!        CODE FOR BOTH INCREMENTS EQUAL TO 1


!        CLEAN-UP LOOP

   20 M = MOD(N,7)
      IF( M == 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DX(I)
   30 CONTINUE
      IF( N < 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        DY(I) = DX(I)
        DY(I + 1) = DX(I + 1)
        DY(I + 2) = DX(I + 2)
        DY(I + 3) = DX(I + 3)
        DY(I + 4) = DX(I + 4)
        DY(I + 5) = DX(I + 5)
        DY(I + 6) = DX(I + 6)
   50 CONTINUE
      RETURN
      END



      REAL ( KIND = wp ) :: FUNCTION DDOT(N,DX,INCX,DY,INCY)

!     FORMS THE DOT PRODUCT OF TWO VECTORS.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.

      INTEGER :: I,INCX,INCY,IX,IY,M,MP1,N
      REAL ( KIND = wp ) :: DX(*),DY(*),DTEMP,ZERO
      PARAMETER ( ZERO = 0.0_wp )

!D    DDOT = ZERO
      DTEMP = ZERO
      IF(N<=0)RETURN
      IF(INCX==1.AND.INCY==1)GO TO 20

!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!          NOT EQUAL TO 1

      IX = 1
      IY = 1
      IF(INCX<0)IX = (-N + 1)*INCX + 1
      IF(INCY<0)IY = (-N + 1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
!D    DDOT = DTEMP
      RETURN

!        CODE FOR BOTH INCREMENTS EQUAL TO 1


!        CLEAN-UP LOOP

   20 M = MOD(N,5)
      IF( M == 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)*DY(I)
   30 CONTINUE
      IF( N < 5 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) + 
     *   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
   50 CONTINUE
   60 CONTINUE
!D    DDOT = DTEMP
      RETURN
      END



      REAL ( KIND = wp ) :: FUNCTION DNRM2 ( N, DX, INCX)
      INTEGER :: N, INCX, NEXT
      REAL ( KIND = wp ) :: DX( * ), CUTLO, CUTHI, HITEST, SUM,
     *                 XMAX, ZERO, ONE
      INTRINSIC        ABS, SQRT
      INTEGER :: I, J, NN
      PARAMETER ( ZERO = 0.0_wp, ONE = 1.0_wp )

!     EUCLIDEAN NORM OF THE N-VECTOR STORED IN DX() WITH STORAGE
!     INCREMENT INCX .
!     IF    N <= 0 RETURN WITH RESULT = 0.
!     IF N >= 1 THEN INCX MUST BE >= 1

!           C.L.LAWSON, 1978 JAN 08

!     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE
!     HOPEFULLY APPLICABLE TO ALL MACHINES.
!         CUTLO = MAXIMUM OF  DSQRT(U/EPS)  OVER ALL KNOWN MACHINES.
!         CUTHI = MINIMUM OF  DSQRT(V)      OVER ALL KNOWN MACHINES.
!     WHERE
!         EPS = SMALLEST NO. SUCH THAT EPS + 1. > 1.
!         U = SMALLEST POSITIVE NO. (UNDERFLOW LIMIT)
!         V = LARGEST  NO. (OVERFLOW  LIMIT)

!     BRIEF OUTLINE OF ALGORITHM..

!     PHASE 1    SCANS ZERO COMPONENTS.
!     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND <= CUTLO
!     MOVE TO PHASE 3 WHEN A COMPONENT IS > CUTLO
!     MOVE TO PHASE 4 WHEN A COMPONENT IS >= CUTHI/M
!     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.

!     VALUES FOR CUTLO AND CUTHI..
!     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER
!     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..
!     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE
!                   UNIVAC AND DEC AT 2**(-103)
!                   THUS CUTLO = 2**(-51) = 4.44089E-16
!     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.
!                   THUS CUTHI = 2**(63.5) = 1.30438E19
!     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC.
!                   THUS CUTLO = 2**(-33.5) = 8.23181D-11
!     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19
!D    DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /

      IF ( N <= 0) THEN
!D       DNRM2 = ZERO
      ELSE
         NEXT = 1
         SUM = ZERO
         NN = N * INCX

!  BEGIN MAIN LOOP

         I = 1
   20    CONTINUE
         GO TO ( 30, 50, 70, 110 ), NEXT
   30    CONTINUE
         IF( ABS( DX( I ) ) > CUTLO ) GO TO 85
         NEXT = 2
         XMAX = ZERO

!  PHASE 1.  SUM IS ZERO

   50    CONTINUE
         IF ( DX( I ) == ZERO ) GO TO 200
         IF ( ABS( DX( I ) ) > CUTLO ) GO TO 85

!  PREPARE FOR PHASE 2.

         NEXT = 3
         GO TO 105

!  PREPARE FOR PHASE 4.

  100    CONTINUE
         I = J
         NEXT = 4
         SUM = ( SUM / DX( I ) ) / DX( I )
  105    CONTINUE
         XMAX = ABS( DX( I ) )
         SUM = SUM + ( DX( I ) / XMAX ) ** 2
         GO TO 200

!  PHASE 2.  SUM IS SMALL. SCALE TO AVOID DESTRUCTIVE UNDERFLOW.

   70    CONTINUE
         IF ( ABS( DX( I ) ) > CUTLO ) THEN

!  PREPARE FOR PHASE 3.

            SUM = ( SUM * XMAX) * XMAX
            GO TO 85
         END IF

!  COMMON CODE FOR PHASES 2 AND 4.
!  IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.

  110    CONTINUE
         IF ( ABS( DX( I ) ) > XMAX ) THEN
            SUM = ONE + SUM * ( XMAX / DX( I ) ) ** 2
            XMAX = ABS( DX( I ) )
         ELSE
            SUM = SUM + ( DX( I ) / XMAX ) ** 2
         END IF
  200    CONTINUE
         I = I + INCX
         IF ( I <= NN ) GO TO 20

!  END OF MAIN LOOP.

!  COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.

!D       DNRM2 = XMAX * SQRT(SUM)
         GO TO 300

!  FOR REAL OR D.P. SET HITEST = CUTHI/N

   85    CONTINUE
         HITEST = CUTHI/FLOAT( N )

!  PHASE 3. SUM IS MID-RANGE.  NO SCALING.

         DO 95 J = I, NN, INCX
            IF( ABS( DX( J ) ) >= HITEST ) GO TO 100
            SUM = SUM + DX( J ) ** 2
   95    CONTINUE
!D       DNRM2 = SQRT( SUM )
      END IF
  300 CONTINUE
      RETURN
      END



      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)

!     CONSTANT TIMES A VECTOR PLUS A VECTOR.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.

      INTEGER :: I,INCX,INCY,IX,IY,M,MP1,N
      REAL ( KIND = wp ) :: DX(*),DY(*),DA,ZERO

!  SET CONSTANT REAL PARAMETERS.

      PARAMETER ( ZERO = 0.0_wp )

      IF(N<=0)RETURN
      IF (DA == ZERO) RETURN
      IF(INCX==1.AND.INCY==1)GO TO 20

!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!          NOT EQUAL TO 1

      IX = 1
      IY = 1
      IF(INCX<0)IX = (-N + 1)*INCX + 1
      IF(INCY<0)IY = (-N + 1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN

!        CODE FOR BOTH INCREMENTS EQUAL TO 1


!        CLEAN-UP LOOP

   20 M = MOD(N,4)
      IF( M == 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF( N < 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
   50 CONTINUE
      RETURN
      END



      SUBROUTINE DROT (N,DX,INCX,DY,INCY,C,S)

!     APPLIES A PLANE ROTATION.
!     JACK DONGARRA, LINPACK, 3/11/78.

      REAL ( KIND = wp ) :: DX(*),DY(*),DTEMP,C,S
      INTEGER :: I,INCX,INCY,IX,IY,N

      IF(N<=0)RETURN
      IF(INCX==1.AND.INCY==1)GO TO 20

!       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
!         TO 1

      IX = 1
      IY = 1
      IF(INCX<0)IX = (-N + 1)*INCX + 1
      IF(INCY<0)IY = (-N + 1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = C*DX(IX) + S*DY(IY)
        DY(IY) = C*DY(IY) - S*DX(IX)
        DX(IX) = DTEMP
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN

!       CODE FOR BOTH INCREMENTS EQUAL TO 1

   20 DO 30 I = 1,N
        DTEMP = C*DX(I) + S*DY(I)
        DY(I) = C*DY(I) - S*DX(I)
        DX(I) = DTEMP
   30 CONTINUE
      RETURN
      END
      SUBROUTINE DROTG(DA,DB,C,S)

!     CONSTRUCT GIVENS PLANE ROTATION.
!     JACK DONGARRA, LINPACK, 3/11/78.

      REAL ( KIND = wp ) :: DA,DB,C,S,ROE,SCALE,R,Z,ZERO,ONE
      INTRINSIC ABS, SQRT, SIGN
      PARAMETER ( ZERO = 0.0_wp, ONE = 1.0_wp )

      ROE = DB
      IF( ABS(DA) > ABS(DB) ) ROE = DA
      SCALE = ABS(DA) + ABS(DB)
      IF( SCALE /= ZERO ) GO TO 10
         C = ONE
         S = ZERO
         R = ZERO
         GO TO 20
   10 R = SCALE*SQRT((DA/SCALE)**2 + (DB/SCALE)**2)
      R = SIGN(ONE,ROE)*R
      C = DA/R
      S = DB/R
   20 Z = ONE
      IF( ABS(DA) > ABS(DB) ) Z = S
      IF( ABS(DB) >= ABS(DA) .AND. C /= ZERO ) Z = ONE/C
      DA = R
      DB = Z
      RETURN
      END



      SUBROUTINE  DSCAL(N,DA,DX,INCX)

!     SCALES A VECTOR BY A CONSTANT.
!     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.

      REAL ( KIND = wp ) :: DA,DX(*)
      INTEGER :: I,INCX,M,MP1,N,NINCX

      IF(N<=0)RETURN
      IF(INCX==1)GO TO 20

!        CODE FOR INCREMENT NOT EQUAL TO 1

      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
        DX(I) = DA*DX(I)
   10 CONTINUE
      RETURN

!        CODE FOR INCREMENT EQUAL TO 1


!        CLEAN-UP LOOP

   20 M = MOD(N,5)
      IF( M == 0 ) GO TO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF( N < 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I + 1) = DA*DX(I + 1)
        DX(I + 2) = DA*DX(I + 2)
        DX(I + 3) = DA*DX(I + 3)
        DX(I + 4) = DA*DX(I + 4)
   50 CONTINUE
      RETURN
      END



      SUBROUTINE  DSWAP (N,DX,INCX,DY,INCY)

!     INTERCHANGES TWO VECTORS.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.

      REAL ( KIND = wp ) :: DX(*),DY(*),DTEMP
      INTEGER :: I,INCX,INCY,IX,IY,M,MP1,N

      IF(N<=0)RETURN
      IF(INCX==1.AND.INCY==1)GO TO 20

!       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
!         TO 1

      IX = 1
      IY = 1
      IF(INCX<0)IX = (-N + 1)*INCX + 1
      IF(INCY<0)IY = (-N + 1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DX(IX)
        DX(IX) = DY(IY)
        DY(IY) = DTEMP
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN

!       CODE FOR BOTH INCREMENTS EQUAL TO 1


!       CLEAN-UP LOOP

   20 M = MOD(N,3)
      IF( M == 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
   30 CONTINUE
      IF( N < 3 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
        DTEMP = DX(I + 1)
        DX(I + 1) = DY(I + 1)
        DY(I + 1) = DTEMP
        DTEMP = DX(I + 2)
        DX(I + 2) = DY(I + 2)
        DY(I + 2) = DTEMP
   50 CONTINUE
      RETURN
      END



!D    INTEGER FUNCTION IDAMAX(N,DX,INCX)

!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)

      REAL ( KIND = wp ) :: DX(*),DMAX
      INTEGER :: I,INCX,IX,N

!D    IDAMAX = 0
      IF( N<1 .OR. INCX<=0 ) RETURN
!D    IDAMAX = 1
      IF(N==1)RETURN
      IF(INCX==1)GO TO 20

!        code for increment not equal to 1

      IX = 1
      DMAX = ABS(DX(1))
      IX = IX + INCX
      DO 10 I = 2,N
         IF(ABS(DX(IX))<=DMAX) GO TO 5
!D       IDAMAX = I
         DMAX = ABS(DX(IX))
    5    IX = IX + INCX
   10 CONTINUE
      RETURN

!        code for increment equal to 1

   20 DMAX = ABS(DX(1))
      DO 30 I = 2,N
         IF(ABS(DX(I))<=DMAX) GO TO 30
!D       IDAMAX = I
         DMAX = ABS(DX(I))
   30 CONTINUE
      RETURN
      END



      REAL ( KIND = wp ) :: FUNCTION DASUM(N,DX,INCX)

!     takes the sum of the absolute values.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)

      REAL ( KIND = wp ) :: DX(*),DTEMP
      INTEGER :: I,INCX,M,MP1,N,NINCX

!D    DASUM = 0.0D0
      DTEMP = 0.0D0
      IF( N<=0 .OR. INCX<=0 )RETURN
      IF(INCX==1)GO TO 20

!        code for increment not equal to 1

      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
        DTEMP = DTEMP + ABS(DX(I))
   10 CONTINUE
!D    DASUM = DTEMP
      RETURN

!        code for increment equal to 1


!        clean-up loop

   20 M = MOD(N,6)
      IF( M == 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + ABS(DX(I))
   30 CONTINUE
      IF( N < 6 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
        DTEMP = DTEMP + ABS(DX(I)) + ABS(DX(I + 1)) + ABS(DX(I + 2))
     * + ABS(DX(I + 3)) + ABS(DX(I + 4)) + ABS(DX(I + 5))
   50 CONTINUE
   60 CONTINUE
!D    DASUM = DTEMP
      RETURN
      END
