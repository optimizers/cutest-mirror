! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE DASMBL( N, NG, MAXSEL, NSEMIB, LH, LIH, NNZH, &
                         NFREE, IFREE, ISTADH, LSTADH, ICNA, LICNA,   &
                         ISTADA, LSTADA, INTVAR, LNTVAR, IELVAR, LELVAR,  &
                         IELING, LELING, ISTADG, LSTADG, ISTAEV, LSTAEV,  &
                         ISTAGV, LNSTGV, ISVGRP, LNVGRP, IRNH, JCNH,    &
                         NXTROW, LNXTRW, IWK, LIWK, A, LA, GUVALS,  &
                         LNGUVL, HUVALS, LNHUVL, GVALS2, GVALS3, GSCALE,  &
                         ESCALE, LESCAL, H, WK, LWK, GXEQX, LGXEQX,  &
                         INTREP, LINTRE, ITYPEE, LITYPE, RANGE, IPRINT,  &
                         IOUT, BAND, MAXSBW, INFORM, NOZERO, FIXSTR)
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  ********************************************************************

!  Assemble the second derivative matrix of a groups partially
!  separable function in either co-ordinate or band format.

!  Nick Gould, February 20TH 1989, for CGT Productions.

!  ********************************************************************

      INTEGER :: N, NG, MAXSEL, NFREE, LNXTRW, NNZH, LITYPE
      INTEGER :: LA, LH, IOUT, INFORM, NSEMIB, IPRINT, LIH
      INTEGER :: LSTADH, LICNA, LSTADA, LNTVAR, LELVAR, LELING
      INTEGER :: LSTADG, LSTAEV, LNSTGV, LNVGRP, LIWK,   MAXSBW
      INTEGER :: LNGUVL, LNHUVL, LESCAL, LWK,    LGXEQX, LINTRE
      LOGICAL :: BAND, NOZERO, FIXSTR
      INTEGER :: IFREE(  N ),      ISTADH( LSTADH ), ICNA( LICNA )  
      INTEGER :: ISTADA( LSTADA ), INTVAR( LNTVAR ), JCNH( LIH )
      INTEGER :: IELVAR( LELVAR ), IELING( LELING ), IRNH( LIH )
      INTEGER :: ISTADG( LSTADG ), ISTAEV( LSTAEV )
      INTEGER :: ISTAGV( LNSTGV ), ISVGRP( LNVGRP ), IWK( LIWK )
      INTEGER :: ITYPEE( LITYPE )
      INTEGER :: NXTROW( 2,        LNXTRW )
      LOGICAL :: GXEQX( LGXEQX ),  INTREP( LINTRE )
      REAL ( KIND = wp ) :: A(     LA ),      GVALS2( NG ),     GVALS3( NG ), &
                       GUVALS( LNGUVL ), HUVALS( LNHUVL ), GSCALE( NG ),  &
                       ESCALE( LESCAL ), H(      LH ),     WK(     LWK )
      EXTERNAL :: RANGE 

!  LOCAL VARIABLES.

      INTEGER :: I, II,  IG, J,  JJ, K,  L, IP,  NN,     NNN
      INTEGER :: NEWPT,  NIN,    IELL,   IEL,    IHNEXT
      INTEGER :: NVAREL, IG1,    LISTVS, LISTVE, IELH
      INTEGER :: INEXT,  IJHESS, IROW,   JCOL,   JCOLST, ISTART
      REAL ( KIND = wp ) :: ONE,    ZERO,   WKI,    HESNEW, GDASH,  G2DASH, &
                       SCALEE
      CHARACTER ( LEN = 2 ) :: MATRIX( 36, 36 )
!D    EXTERNAL         DSETVL, DSETVI
      INTRINSIC        ABS,    MAX,    MIN

!  Set constant real parameters.

      PARAMETER ( ZERO = 0.0_wp, ONE = 1.0_wp )
      IF ( .NOT. BAND .AND. NFREE > LNXTRW ) GO TO 620

!  Renumber the free variables so that they are variables 1 to NFREE.

      DO 10 I = 1, N
         IWK( I ) = 0
   10 CONTINUE
      DO 20 I = 1, NFREE
         IWK( IFREE( I ) ) = I

!  Initialize the link list which points to the row numbers which
!  are used in the columns of the assembled Hessian.

!  NXTROW( 1, . ) gives the link list. the list for column J starts
!                 in NXTROW( 1, J ) and ends when NXTROW( 1, K ) = - 1.
!  NXTROW( 2, . ) gives the position in H of the current link.

         IF ( .NOT. BAND ) NXTROW( 1, I ) = - 1
   20 CONTINUE
      IF ( IPRINT >= 10 ) THEN
         WRITE( IOUT, 2060 ) NFREE, ( IFREE( I ), I = 1, NFREE )
      END IF

!  If a band storage scheme is to be used, initialize the entries
!  within the band as zero.

      IF ( BAND ) THEN
         MAXSBW = 0
         IF ( NFREE * ( NSEMIB + 1 ) > LH ) GO TO 610
         DO 30 I = 1, NFREE * ( NSEMIB + 1 )
            H( I ) = ZERO
   30    CONTINUE   
      ELSE
         NNZH = 0
         NEWPT = NFREE + 1
         IF ( NEWPT > LNXTRW ) GO TO 620
      END IF

! ------------------------------------------------------
!  Form the rank-one second order term for the Ith group.
! ------------------------------------------------------

      DO 200 IG = 1, NG
         IF ( IPRINT >= 100 ) WRITE( IOUT, 2070 ) IG
         IF ( GXEQX( IG ) ) GO TO 200
         IF ( .NOT. FIXSTR .AND. GSCALE( IG ) == ZERO ) GO TO 200
         G2DASH = GSCALE( IG ) * GVALS3( IG )
         IF ( IPRINT >= 100 ) WRITE( 6, * ) ' GVALS3(IG) ', GVALS3(IG)
         IF ( NOZERO .AND. G2DASH == ZERO ) GO TO 200
         IG1 = IG + 1
         LISTVS = ISTAGV( IG )
         LISTVE = ISTAGV( IG1 ) - 1

!  Form the gradient of the IG-th group.

      CALL DSETVI( LISTVE - LISTVS + 1, WK, ISVGRP( LISTVS ), ZERO )

!  Consider any nonlinear elements for the group.

         DO 130 IELL = ISTADG( IG ), ISTADG( IG1 ) - 1
            IEL = IELING( IELL )
            K = INTVAR( IEL )
            L = ISTAEV( IEL )
            NVAREL = ISTAEV( IEL + 1 ) - L
            SCALEE = ESCALE( IELL )
            IF ( INTREP( IEL ) ) THEN

!  The IEL-th element has an internal representation.

               NIN = INTVAR( IEL + 1 ) - K
               CALL RANGE ( IEL, .TRUE., GUVALS( K ), &
                            WK( N + 1 ), NVAREL, NIN,  &
                            ITYPEE( IEL ), NIN, NVAREL )
               DO 110 I = 1, NVAREL
                  J = IELVAR( L )
                  WK( J ) = WK( J ) + SCALEE * WK( N + I )
                  L = L + 1
  110          CONTINUE
            ELSE

!  The IEL-th element has no internal representation.

               DO 120 I = 1, NVAREL
                  J = IELVAR( L )
                  WK( J ) = WK( J ) + SCALEE * GUVALS( K )
                  K = K + 1
                  L = L + 1
  120          CONTINUE
            END IF
  130    CONTINUE

!  Include the contribution from the linear element.

         DO 140 K = ISTADA( IG ), ISTADA( IG1 ) - 1
            J = ICNA( K )
            WK( J ) = WK( J ) + A( K )
  140    CONTINUE

!  The gradient is complete. form the J-th column of the rank-one matrix

         DO 190 L = LISTVS, LISTVE
            JJ = ISVGRP( L )
            J = IWK( JJ )
            IF ( J == 0 ) GO TO 190

!  Find the entry in row i of this column.

            DO 180 K = LISTVS, LISTVE
               II = ISVGRP( K )
               I = IWK( II )
               IF ( I == 0 .OR. I > J ) GO TO 180

!  Skip all elements which lie outside a band of width NSEMIB.

               IF ( BAND ) MAXSBW = MAX( MAXSBW, J - I )
               IF ( J - I > NSEMIB ) GO TO 180
               HESNEW = WK( II ) * WK( JJ ) * G2DASH
               IF ( IPRINT >= 100 ) WRITE( IOUT, 2090 ) I, J, HESNEW
               IF ( NOZERO .AND. HESNEW == ZERO ) GO TO 180

!  Obtain the appropriate storage location in H for the new entry.


!  Case 1: band matrix storage scheme.

               IF ( BAND ) THEN

!  The entry belongs on the diagonal.

                  IF ( I == J ) THEN
                     H( I ) = H( I ) + HESNEW

!  The entry belongs off the diagonal.

                  ELSE
                     H( NFREE + J - I + NSEMIB * ( I - 1 ) ) =  &
               H( NFREE + J - I + NSEMIB * ( I - 1 ) ) + HESNEW
                  END IF   

!  Case 2: coordinate storage scheme.

               ELSE
                  ISTART = J
  150             CONTINUE
                  INEXT = NXTROW( 1, ISTART )
                  IF ( INEXT == - 1 ) THEN
                     IF ( NEWPT > LNXTRW ) GO TO 620

!  The (I,J)-th location is empty. Place the new entry in this location
!  and add another link to the list.

                     NNZH = NNZH + 1
                     IF ( NNZH > LH .OR. NNZH > LIH ) GO TO 610
                     IRNH( NNZH ) = I
                     JCNH( NNZH ) = J
                     H( NNZH ) = HESNEW
                     NXTROW( 1, ISTART ) = NEWPT
                     NXTROW( 2, ISTART ) = NNZH
                     NXTROW( 1, NEWPT ) = - 1
                     NEWPT = NEWPT + 1
                  ELSE

!  Continue searching the linked list for an entry in row I, column J.

                     IF ( IRNH( NXTROW( 2, ISTART ) ) == I ) THEN
                        IP = NXTROW( 2, ISTART )
                        H( IP ) = H( IP ) + HESNEW
                     ELSE
                        ISTART = INEXT
                        GO TO 150
                     END IF
                  END IF
               END IF
  180       CONTINUE
  190    CONTINUE
  200 CONTINUE

!  Reset the workspace array to zero.

      CALL DSETVL( MAXSEL, WK, 1, ZERO )

! --------------------------------------------------------
!  Add on the low rank first order terms for the I-th group.
! --------------------------------------------------------

      DO 300 IG = 1, NG
         IF ( IPRINT >= 100 ) WRITE( IOUT, 2100 ) IG
         IF ( .NOT. FIXSTR .AND. GSCALE( IG ) == ZERO ) GO TO 300
         IF ( GXEQX( IG ) ) THEN
            GDASH = GSCALE( IG )
         ELSE
            GDASH = GSCALE( IG ) * GVALS2( IG )
            IF ( NOZERO .AND. GDASH == ZERO ) GO TO 300
         END IF
         IG1 = IG + 1

!  See if the group has any nonlinear elements.

         DO 290 IELL = ISTADG( IG ), ISTADG( IG1 ) - 1
            IEL = IELING( IELL )
            LISTVS = ISTAEV( IEL )
            LISTVE = ISTAEV( IEL + 1 ) - 1
            NVAREL = LISTVE - LISTVS + 1
            IELH = ISTADH( IEL )
            IHNEXT = IELH
            SCALEE = ESCALE( IELL )
            DO 250 L = LISTVS, LISTVE
               J = IWK( IELVAR( L ) )
               IF ( J /= 0 ) THEN

!  The IEL-th element has an internal representation.
!  Compute the J-th column of the element Hessian matrix.

                  IF ( INTREP( IEL ) ) THEN

!  Compute the J-th column of the Hessian.

                     WK( L - LISTVS + 1 ) = ONE

!  Find the internal variables.

                     NIN = INTVAR( IEL + 1 ) - INTVAR( IEL )
                     CALL RANGE ( IEL, .FALSE., WK( 1 ), &
                              WK( MAXSEL + 1 ), NVAREL, NIN, &
                              ITYPEE( IEL ), NVAREL, NIN )

!  Multiply the internal variables by the element Hessian.

                     NN = MAXSEL + NIN
      CALL DSETVL( NIN, WK( NN + 1 ), 1, ZERO )

!  Only the upper triangle of the element Hessian is stored.

                     JCOLST = IELH - 1
                     DO 220 JCOL = 1, NIN
                        IJHESS = JCOLST
                        JCOLST = JCOLST + JCOL
                        WKI = WK( MAXSEL + JCOL ) * GDASH
                        DO 210 IROW = 1, NIN
                           IF ( IROW <= JCOL ) THEN
                              IJHESS = IJHESS + 1
                           ELSE
                              IJHESS = IJHESS + IROW - 1
                           END IF
                           WK( NN + IROW ) = WK( NN + IROW ) + &
                                             WKI * HUVALS( IJHESS )
  210                   CONTINUE
  220                CONTINUE

!  Scatter the product back onto the elemental variables.

                     NNN = NN + NIN
                     CALL RANGE ( IEL, .TRUE., WK( NN + 1 ), &
                                  WK( NNN + 1 ), NVAREL, NIN, &
                                  ITYPEE( IEL ), NIN, NVAREL )
                     WK( L - LISTVS + 1 ) = ZERO
                  END IF

!  Find the entry in row I of this column.

                  DO 240 K = LISTVS, L
                     I = IWK( IELVAR( K ) )

!  Skip all elements which lie outside a band of width NSEMIB.

                     IF ( BAND .AND. I /= 0 ) MAXSBW =  &
                          MAX( MAXSBW, ABS( J - I ) )
                     IF ( ABS( I - J ) <= NSEMIB .AND. I /= 0 ) THEN

!  Only the upper triangle of the matrix is stored.

                        IF ( I <= J ) THEN
                           II = I
                           JJ = J
                        ELSE
                           II = J
                           JJ = I
                        END IF

!  Obtain the appropriate storage location in H for the new entry.

                        IF ( INTREP( IEL ) ) THEN
                           HESNEW = SCALEE * WK( NNN + K - LISTVS + 1 )
                        ELSE
                           HESNEW = SCALEE * HUVALS( IHNEXT ) * GDASH
                        END IF
                        IF ( IPRINT >= 100 ) &
                           WRITE( 6, 2080 ) II, JJ, IEL, HESNEW

!  Case 1: band matrix storage scheme.

                        IF ( BAND ) THEN

!  The entry belongs on the diagonal.

                           IF ( II == JJ ) THEN
                              H( II ) = H( II ) + HESNEW
                              IF ( K /= L ) H( II ) = H( II ) + HESNEW

!  The entry belongs off the diagonal.

                           ELSE
                              H( NFREE + JJ - II + NSEMIB * ( II - 1 ) ) &
 =  H( NFREE + JJ - II + NSEMIB * ( II - 1 ) ) &
 + HESNEW
                           END IF   

!  Case 2: coordinate storage scheme.

                        ELSE
                           IF ( .NOT. NOZERO .OR. HESNEW /= ZERO)THEN
                              ISTART = JJ
  230                         CONTINUE
                              INEXT = NXTROW( 1, ISTART )
                              IF ( INEXT == - 1 ) THEN
                                 IF ( NEWPT > LNXTRW ) GO TO 620

!  The (I,J)-th location is empty. Place the new entry in this location
!  and add another link to the list.

                                 NNZH = NNZH + 1
                                 IF ( NNZH > LH .OR.  &
                                      NNZH > LIH ) GO TO 610
                                 IRNH( NNZH ) = II
                                 JCNH( NNZH ) = JJ
                                 H( NNZH ) = HESNEW
                                 IF ( K /= L .AND. II == JJ )   &
                                    H( NNZH ) = HESNEW + HESNEW
                                 NXTROW( 1, ISTART ) = NEWPT
                                 NXTROW( 2, ISTART ) = NNZH
                                 NXTROW( 1, NEWPT ) = - 1
                                 NEWPT = NEWPT + 1
                              ELSE

!  Continue searching the linked list for an entry in row I, column J.

                              IF ( IRNH( NXTROW( 2, ISTART ) ) &
 == II ) THEN
                                    IP = NXTROW( 2, ISTART )
                                    H( IP ) = H( IP ) + HESNEW
                                    IF ( K /= L .AND. II == JJ )   &
                                       H( IP ) = H( IP ) + HESNEW
                                 ELSE
                                    ISTART = INEXT
                                    GO TO 230
                                 END IF
                              END IF
                           END IF
                        END IF
                     END IF
                     IHNEXT = IHNEXT + 1
  240             CONTINUE
               ELSE
                  IHNEXT = IHNEXT + L - LISTVS + 1
               END IF
  250       CONTINUE
  290    CONTINUE
  300 CONTINUE

! ----------------------------------------
!  For debugging, print the nonzero values.
! ----------------------------------------

      IF ( IPRINT >= 10 ) THEN
         IF ( .NOT. BAND ) WRITE( IOUT, 2000 ) &
 ( IRNH( I ), JCNH( I ), H( I ), I = 1, NNZH )
         IF ( NFREE <= 36 ) THEN
            DO 420 I = 1, NFREE
               DO 410 J = 1, NFREE
                  MATRIX( I, J ) = '  '
  410          CONTINUE
  420       CONTINUE
            IF ( BAND ) THEN
               DO 440 I = 1, NFREE
                  IF ( H( I ) /= ZERO ) MATRIX( I, I ) = ' *'
                  DO 430 J = 1, MIN( NSEMIB, NFREE - I )
                     K = NFREE + J + NSEMIB * ( I - 1 )
                     IF ( H( K ) /= ZERO ) THEN
                        MATRIX( I + J, I ) = ' *'
                        MATRIX( I, I + J ) = ' *'
                     END IF
  430             CONTINUE
  440          CONTINUE
            ELSE
               DO 450 I = 1, NNZH
                  IF ( IRNH( I ) > NFREE ) THEN
                     WRITE( IOUT, * ) ' ENTRY OUT OF BOUNDS IN ASEMBL ', &
                                      ' ROW NUMBER = ', IRNH( I )
                     STOP
                  END IF
                  IF ( JCNH( I ) > NFREE ) THEN
                     WRITE( IOUT, * ) ' ENTRY OUT OF BOUNDS IN ASEMBL ', &
                                      ' COL NUMBER = ', JCNH( I )
                     STOP
                  END IF
                  MATRIX( IRNH( I ), JCNH( I ) ) = ' *'
                  MATRIX( JCNH( I ), IRNH( I ) ) = ' *'
  450          CONTINUE
            END IF
            WRITE( IOUT, 2040 ) ( I, I = 1, NFREE )
            DO 460 I = 1, NFREE
               WRITE( IOUT, 2050 ) I, ( MATRIX( I, J ), J = 1, NFREE )
  460       CONTINUE
         END IF
      END IF
      INFORM = 0
      RETURN

!  Unsuccessful returns.

  610 CONTINUE
      INFORM = 1
      IF ( LH <= LIH ) THEN
         WRITE( IOUT, 2010 ) LH
      ELSE
         WRITE( IOUT, 2030 ) LIH
      END IF
      RETURN
  620 CONTINUE
      INFORM = 2
      WRITE( IOUT, 2020 ) LNXTRW
      RETURN

!  Non-executable statements

 2000 FORMAT( '    Row  Column    Value        Row  Column    Value ', / &
              ' --- ------ ----- --- ------ ----- ', / &
 ( 2I6, 1P, D24.16, 2I6, 1P, D24.16 ) )
 2010 FORMAT( ' ** Array dimension LH =', I6, ' too small in ASMBL. ' )
 2020 FORMAT( ' ** Array dimension LNXTRW =',I8, ' too small in ASMBL.')
 2030 FORMAT( ' ** Array dimension LIH =', I6, ' too small in ASMBL.' )
 2040 FORMAT( /, 5X, 36I2 )
 2050 FORMAT( I3, 2X, 36A2 )
 2060 FORMAT( /, I6, ' free variables. They are ', 8I5, /, ( 14I5 ) )
 2070 FORMAT( ' Group ', I5, ' rank-one terms ' )
 2080 FORMAT( ' Row ', I6, ' Column ', I6, ' used from element ', I6, &
              ' value = ', 1P, D24.16 )
 2090 FORMAT( ' Row ', I6, ' column ', I6, ' used. Value = ',1P,D24.16)
 2100 FORMAT( ' Group ', I5, ' second-order terms ' )

!  END OF ASMBL.

      END
