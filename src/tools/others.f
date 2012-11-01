! ( Last modified on 23 Dec 2000 at 22:01:38 )



      SUBROUTINE DGELIM( M, N, IPVT, JCOL, A )
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: M, N
      INTEGER :: IPVT ( M ), JCOL ( N )
      REAL ( KIND = wp ) :: A ( M, N )

!  Perform the first M steps of Gaussian Elimination with
!  complete pivoting on the M by N ( M <= N) matrix A.

!  Nick Gould, 23rd September 1991.
!  For CGT productions.

      INTEGER :: I, J, K, IPIVOT, JPIVOT
      REAL ( KIND = wp ) :: APIVOT, ONE, ATEMP
      PARAMETER ( ONE = 1.0_wp )

!  Initialize the column indices.

      DO 10 J = 1, N
         JCOL( J ) = J
   10 CONTINUE

!  Main loop.

      DO 100 K = 1, M

!  Compute the K-th pivot.

         APIVOT = - ONE
         DO 30 J = K, N
            DO 20 I = K, M
               IF ( ABS( A( I, J ) ) > APIVOT ) THEN
                  APIVOT = ABS( A( I, J ) )
                  IPIVOT = I
                  JPIVOT = J
               END IF
   20       CONTINUE
   30    CONTINUE

!  Interchange rows I and IPIVOT.

         IPVT( K ) = IPIVOT
         IF ( IPIVOT > K ) THEN
            DO 40 J = K, N
               ATEMP = A( IPIVOT, J )
               A( IPIVOT, J ) = A( K, J )
               A( K, J ) = ATEMP
   40       CONTINUE
         END IF

!  Interchange columns J and JPIVOT.

         IF ( JPIVOT > K ) THEN
            J = JCOL( JPIVOT )
            JCOL( JPIVOT ) = JCOL( K )
            JCOL( K ) = J
            DO 50 I = 1, M
               ATEMP = A( I, JPIVOT )
               A( I, JPIVOT ) = A( I, K )
               A( I, K ) = ATEMP
   50       CONTINUE
         END IF

!  Perform the elimination.

         APIVOT = A( K, K )
         DO 70 I = K + 1, M
            ATEMP = A( I, K ) / APIVOT
            A( I, K ) = ATEMP
            DO 60 J = K + 1, N
               A( I, J ) = A( I, J ) - ATEMP * A( K, J )
   60       CONTINUE
   70    CONTINUE
  100 CONTINUE
      RETURN

!  End of subroutine GELIM.

      END



      SUBROUTINE DGESLV( M, IPVT, A, X )
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: M
      INTEGER :: IPVT ( M )
      REAL ( KIND = wp ) :: A ( M, M ), X ( M )

!  Solve the equations A(T)x = b. The vector b is input in X.
!  The LU factors of P A are input in A; The permutation P is stored
!  in IPVT. The solution x is output in X.

!  Nick Gould, 23rd September 1991.
!  For CGT productions.

      INTEGER :: I, K
      REAL ( KIND = wp ) :: XTEMP, ZERO
      PARAMETER ( ZERO = 0.0_wp )

!  Solve U(T)y = b. The vector b is input in X; y is output in X.

      DO 20 K = 1, M
         XTEMP = ZERO
         DO 10 I = 1, K - 1
            XTEMP = XTEMP + A( I, K ) * X( I )
   10    CONTINUE
         X( K ) = ( X( K ) - XTEMP ) / A( K, K )
   20 CONTINUE

!  Solve L(T) x = y. The vector y is input in X; x is output in X.

      DO 40 K = M - 1, 1, - 1
         XTEMP = ZERO
         DO 30 I = K + 1, M
            XTEMP = XTEMP + A( I, K ) * X( I )
   30    CONTINUE
         X( K ) = X( K ) - XTEMP
         I = IPVT( K )
         IF ( I /= K ) THEN
            XTEMP = X( I )
            X( I ) = X( K )
            X( K ) = XTEMP
         END IF
   40 CONTINUE
      RETURN

!  End of subroutine GESLV.

      END


!  ** FOR THE CRAY 2, LINES STARTING 'CDIR$ IVDEP' TELL THE COMPILER TO
!     IGNORE POTENTIAL VECTOR DEPENDENCIES AS THEY ARE KNOWN TO BE O.K.

      SUBROUTINE DSYMMH( MAXSZH, ISYMMH, ISYMMD )
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: MAXSZH
      INTEGER :: ISYMMH( MAXSZH, MAXSZH ), ISYMMD( MAXSZH )

!  GIVEN A COLUMNWISE STORAGE SCHEME OF THE UPPER TRIANGLE OF A
!  SYMMETRIC MATRIX OF ORDER MAXSZH, COMPUTE THE POSITION OF THE
!  I,J-TH ENTRY OF THE SYMMETRIC MATRIX IN THIS SCHEME.

!  THE VALUE ISYMMH( I, J ) + 1 GIVES THE POSITION OF THE I,J-TH
!  ENTRY OF THE MATRIX IN THE UPPER TRIANGULAR SCHEME.

!  NICK GOULD, 10TH OF MAY 1989.
!  FOR CGT PRODUCTIONS.

      INTEGER :: I, J, K
      K = 0
      DO 20 J = 1, MAXSZH
!DIR$ IVDEP
         DO 10 I = 1, J - 1
            ISYMMH( I, J ) = K
            ISYMMH( J, I ) = K
            K = K + 1
   10    CONTINUE
         ISYMMD( J ) = K
         ISYMMH( J, J ) = K
         K = K + 1
   20 CONTINUE
      RETURN

!  END OF SYMMH.

      END



      SUBROUTINE DSETVL( N, X, INCX, VL )

!     ******************************************************************

      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: INCX, N
      REAL ( KIND = wp ) :: VL, X( * )
      INTEGER :: IX, M, I, J
      IF ( N <= 0 ) RETURN
      IF ( INCX /= 1 ) THEN
         IF( INCX < 0 ) THEN
             IX = ( - N + 1 ) * INCX + 1
         ELSE
             IX = 1
         END IF
         J = IX + ( N - 1 ) * INCX
         DO 100 I = IX, J, INCX
            X( I ) = VL
  100    CONTINUE
      ELSE
         M = MOD( N, 5 )
         IF ( M /= 0 ) THEN
            DO 200 I = 1, M
               X( I ) = VL
  200       CONTINUE
         END IF
         IF ( N >= 5 ) THEN
            IX = M + 1
            DO 300 I = IX, N, 5
               X( I ) = VL
               X( I + 1 ) = VL
               X( I + 2 ) = VL
               X( I + 3 ) = VL
               X( I + 4 ) = VL
  300       CONTINUE
         END IF
      END IF
      RETURN
      END



      SUBROUTINE DSETVI( NVAR, X, IVAR, VL )

!     ******************************************************************

      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: NVAR, IVAR( * ), I
      REAL ( KIND = wp ) :: VL, X( * )
      IF ( NVAR <= 0 ) RETURN
      DO 10 I = 1, NVAR
         X( IVAR( I ) ) = VL
   10 CONTINUE
      RETURN
      END
