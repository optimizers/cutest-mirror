! ( Last modified on 23 Dec 2000 at 22:01:38 )

      SUBROUTINE GELIM( m, n, IPVT, jcol, A )
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: m, n
      INTEGER :: IPVT ( m ), jcol ( n )
      REAL ( KIND = wp ) :: A ( m, n )

!  Perform the first m steps of Gaussian Elimination with
!  complete pivoting on the m by n ( m <= N) matrix A.

!  Nick Gould, 23rd September 1991.
!  For CGT productions.

      INTEGER :: i, j, k, ipivot, jpivot
      REAL ( KIND = wp ) :: apivot, one, atemp
      PARAMETER ( one = 1.0_wp )

!  Initialize the column indices.

      DO 10 j = 1, n
         JCOL( j ) = j
   10 CONTINUE

!  Main loop.

      DO 100 k = 1, m

!  Compute the K-th pivot.

         apivot = - one
         DO 30 j = k, n
            DO 20 i = k, m
               IF ( ABS( A( i, j ) ) > apivot ) THEN
                  apivot = ABS( A( i, j ) )
                  ipivot = i
                  jpivot = j
               END IF
   20       CONTINUE
   30    CONTINUE

!  Interchange rows i and IPIVOT.

         IPVT( k ) = ipivot
         IF ( ipivot > k ) THEN
            DO 40 j = k, n
               atemp = A( ipivot, j )
               A( ipivot, j ) = A( k, j )
               A( k, j ) = atemp
   40       CONTINUE
         END IF

!  Interchange columns j and JPIVOT.

         IF ( jpivot > k ) THEN
            j = JCOL( jpivot )
            JCOL( jpivot ) = JCOL( k )
            JCOL( k ) = j
            DO 50 i = 1, m
               atemp = A( i, jpivot )
               A( i, jpivot ) = A( i, k )
               A( i, k ) = atemp
   50       CONTINUE
         END IF

!  Perform the elimination.

         apivot = A( k, k )
         DO 70 i = k + 1, m
            atemp = A( i, k ) / apivot
            A( i, k ) = atemp
            DO 60 j = k + 1, n
               A( i, j ) = A( i, j ) - atemp * A( k, j )
   60       CONTINUE
   70    CONTINUE
  100 CONTINUE
      RETURN

!  End of subroutine GELIM.

      END SUBROUTINE GELIM

      SUBROUTINE GESLV( m, IPVT, A, X )
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: m
      INTEGER :: IPVT ( m )
      REAL ( KIND = wp ) :: A ( m, m ), X ( m )

!  Solve the equations A(T)x = b. The vector b is input in X.
!  The lu factors of P A are input in A; The permutation P is stored
!  in IPVT. The solution x is output in X.

!  Nick Gould, 23rd September 1991.
!  For CGT productions.

      INTEGER :: i, k
      REAL ( KIND = wp ) :: xtemp, zero
      PARAMETER ( zero = 0.0_wp )

!  Solve U(T)y = b. The vector b is input in X; y is output in X.

      DO 20 k = 1, m
         xtemp = zero
         DO 10 i = 1, k - 1
            xtemp = xtemp + A( i, k ) * X( i )
   10    CONTINUE
         X( k ) = ( X( k ) - xtemp ) / A( k, k )
   20 CONTINUE

!  Solve L(T) x = y. The vector y is input in X; x is output in X.

      DO 40 k = m - 1, 1, - 1
         xtemp = zero
         DO 30 i = k + 1, m
            xtemp = xtemp + A( i, k ) * X( i )
   30    CONTINUE
         X( k ) = X( k ) - xtemp
         i = IPVT( k )
         IF ( i /= k ) THEN
            xtemp = X( i )
            X( i ) = X( k )
            X( k ) = xtemp
         END IF
   40 CONTINUE
      RETURN

!  End of subroutine GESLV.

      END SUBROUTINE GESLV


!  ** FOR THE CRAY 2, LINES STARTING 'CDIR$ IVDEP' TELL THE COMPILER TO
!     IGNORE POTENTIAL VECTOR DEPENDENCIES AS THEY ARE KNOWN TO BE O.K.

      SUBROUTINE SYMMH( maxszh, ISYMMH, ISYMMD )
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: maxszh
      INTEGER :: ISYMMH( maxszh, maxszh ), ISYMMD( maxszh )

!  GIVEN A COLUMNWISE STORAGE SCHEME OF THE UPPER TRIANGLE OF A
!  SYMMETRIC MATRIX OF ORDER maxszh, COMPUTE THE POSITION OF THE
!  i,J-TH ENTRY OF THE SYMMETRIC MATRIX IN THIS SCHEME.

!  THE VALUE ISYMMH( i, j ) + 1 GIVES THE POSITION OF THE i,J-TH
!  ENTRY OF THE MATRIX IN THE UPPER TRIANGULAR SCHEME.

!  NICK GOULD, 10TH OF MAY 1989.
!  FOR CGT PRODUCTIONS.

      INTEGER :: i, j, k
      k = 0
      DO 20 j = 1, maxszh
!DIR$ IVDEP
         DO 10 i = 1, j - 1
            ISYMMH( i, j ) = k
            ISYMMH( j, i ) = k
            k = k + 1
   10    CONTINUE
         ISYMMD( j ) = k
         ISYMMH( j, j ) = k
         k = k + 1
   20 CONTINUE
      RETURN

!  END OF SYMMH.

      END SUBROUTINE SYMMH



      SUBROUTINE SETVL( n, X, incx, VL )

!     ******************************************************************

      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: incx, n
      REAL ( KIND = wp ) :: VL, X( * )
      INTEGER :: ix, m, i, j
      IF ( n <= 0 ) RETURN
      IF ( incx /= 1 ) THEN
         IF( incx < 0 ) THEN
             ix = ( - n + 1 ) * incx + 1
         ELSE
             ix = 1
         END IF
         j = ix + ( n - 1 ) * incx
         DO 100 i = ix, j, incx
            X( i ) = VL
  100    CONTINUE
      ELSE
         m = MOD( n, 5 )
         IF ( m /= 0 ) THEN
            DO 200 i = 1, m
               X( i ) = VL
  200       CONTINUE
         END IF
         IF ( n >= 5 ) THEN
            ix = m + 1
            DO 300 i = ix, n, 5
               X( i ) = VL
               X( i + 1 ) = VL
               X( i + 2 ) = VL
               X( i + 3 ) = VL
               X( i + 4 ) = VL
  300       CONTINUE
         END IF
      END IF
      RETURN
      END SUBROUTINE SETVL



      SUBROUTINE SETVI( nvar, X, IVAR, VL )

!     ******************************************************************

      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: nvar, IVAR( * ), i
      REAL ( KIND = wp ) :: VL, X( * )
      IF ( nvar <= 0 ) RETURN
      DO 10 i = 1, nvar
         X( IVAR( i ) ) = VL
   10 CONTINUE
      RETURN
      END SUBROUTINE SETVI
