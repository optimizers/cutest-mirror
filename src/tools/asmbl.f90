! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE ASMBL( n, ng, maxsel, nsemib, lh, lih, nnzh, &
                         nfree, IFREE, ISTADH, lstadh, ICNA, licna,   &
                         ISTADA, lstada, INTVAR, lntvar, IELVAR, lelvar,  &
                         IELING, leling, ISTADG, lstadg, ISTAEV, lstaev,  &
                         ISTAGV, lnstgv, ISVGRP, lnvgrp, irnh, JCNH,    &
                         NXTROW, lnxtrw, IWK, liwk, A, la, GUVALS,  &
                         lnguvl, HUVALS, lnhuvl, GVALS2, GVALS3, GSCALE,  &
                         ESCALE, lescal, H, WK, lwk, gxeqx, lgxeqx,  &
                         intrep, lintre, ITYPEE, litype, RANGE, iprint,  &
                         iout, BAND, maxsbw, inform, NOZERO, FIXSTR)
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  ********************************************************************

!  Assemble the second derivative matrix of a groups partially
!  separable function in either co-ordinate or band format.

!  Nick Gould, February 20TH 1989, for CGT Productions.

!  ********************************************************************

      INTEGER :: n, ng, maxsel, nfree, lnxtrw, nnzh, litype
      INTEGER :: la, lh, iout, inform, nsemib, iprint, lih
      INTEGER :: lstadh, licna, lstada, lntvar, lelvar, leling
      INTEGER :: lstadg, lstaev, lnstgv, lnvgrp, liwk,   maxsbw
      INTEGER :: lnguvl, lnhuvl, lescal, lwk,    lgxeqx, lintre
      LOGICAL :: BAND, NOZERO, FIXSTR
      INTEGER :: IFREE(  n ),      ISTADH( lstadh ), ICNA( licna )  
      INTEGER :: ISTADA( lstada ), INTVAR( lntvar ), JCNH( lih )
      INTEGER :: IELVAR( lelvar ), IELING( leling ), IRNH( lih )
      INTEGER :: ISTADG( lstadg ), ISTAEV( lstaev )
      INTEGER :: ISTAGV( lnstgv ), ISVGRP( lnvgrp ), IWK( liwk )
      INTEGER :: ITYPEE( litype )
      INTEGER :: NXTROW( 2,        lnxtrw )
      LOGICAL :: GXEQX( lgxeqx ),  INTREP( lintre )
      REAL ( KIND = wp ) :: A(     la ),      GVALS2( ng ),     GVALS3( ng ), &
                       GUVALS( lnguvl ), HUVALS( lnhuvl ), GSCALE( ng ),  &
                       ESCALE( lescal ), H(      lh ),     WK(     lwk )
      EXTERNAL :: RANGE 

!  LOCAL VARIABLES.

      INTEGER :: i, ii,  ig, j,  jj, k,  l, ip,  nn,     nnn
      INTEGER :: newpt,  nin,    iell,   iel,    ihnext
      INTEGER :: nvarel, ig1,    listvs, listve, ielh
      INTEGER :: inext,  ijhess, irow,   jcol,   jcolst, istart
      REAL ( KIND = wp ) :: one,    zero,   wki,    hesnew, gdash,  g2dash, &
                       scalee
      CHARACTER ( LEN = 2 ) :: MATRIX( 36, 36 )
!D    EXTERNAL         SETVL, SETVI
      INTRINSIC        ABS,    MAX,    MIN

!  Set constant real parameters.

      PARAMETER ( zero = 0.0_wp, one = 1.0_wp )
      IF ( .NOT. BAND .AND. nfree > lnxtrw ) GO TO 620

!  Renumber the free variables so that they are variables 1 to NFREE.

      DO 10 i = 1, n
         IWK( i ) = 0
   10 CONTINUE
      DO 20 i = 1, nfree
         IWK( IFREE( i ) ) = i

!  Initialize the link list which points to the row numbers which
!  are used in the columns of the assembled Hessian.

!  NXTROW( 1, . ) gives the link list. the list for column j starts
!                 in NXTROW( 1, j ) and ends when NXTROW( 1, k ) = - 1.
!  NXTROW( 2, . ) gives the position in H of the current link.

         IF ( .NOT. BAND ) NXTROW( 1, i ) = - 1
   20 CONTINUE
      IF ( iprint >= 10 ) THEN
         WRITE( iout, 2060 ) nfree, ( IFREE( i ), i = 1, nfree )
      END IF

!  If a band storage scheme is to be used, initialize the entries
!  within the band as zero.

      IF ( BAND ) THEN
         maxsbw = 0
         IF ( nfree * ( nsemib + 1 ) > lh ) GO TO 610
         DO 30 i = 1, nfree * ( nsemib + 1 )
            H( i ) = zero
   30    CONTINUE   
      ELSE
         nnzh = 0
         newpt = nfree + 1
         IF ( newpt > lnxtrw ) GO TO 620
      END IF

! ------------------------------------------------------
!  Form the rank-one second order term for the Ith group.
! ------------------------------------------------------

      DO 200 ig = 1, ng
         IF ( iprint >= 100 ) WRITE( iout, 2070 ) ig
         IF ( GXEQX( ig ) ) GO TO 200
         IF ( .NOT. FIXSTR .AND. GSCALE( ig ) == zero ) GO TO 200
         g2dash = GSCALE( ig ) * GVALS3( ig )
         IF ( iprint >= 100 ) WRITE( 6, * ) ' GVALS3(IG) ', GVALS3(IG)
         IF ( NOZERO .AND. g2dash == zero ) GO TO 200
         ig1 = ig + 1
         listvs = ISTAGV( ig )
         listve = ISTAGV( ig1 ) - 1

!  Form the gradient of the IG-th group.

      CALL SETVI( listve - listvs + 1, WK, ISVGRP( listvs ), zero )

!  Consider any nonlinear elements for the group.

         DO 130 iell = ISTADG( ig ), ISTADG( ig1 ) - 1
            iel = IELING( iell )
            k = INTVAR( iel )
            l = ISTAEV( iel )
            nvarel = ISTAEV( iel + 1 ) - l
            scalee = ESCALE( iell )
            IF ( INTREP( iel ) ) THEN

!  The IEL-th element has an internal representation.

               nin = INTVAR( iel + 1 ) - k
               CALL RANGE ( iel, .TRUE., GUVALS( k ), &
                            WK( n + 1 ), nvarel, nin,  &
                            ITYPEE( iel ), nin, nvarel )
               DO 110 i = 1, nvarel
                  j = IELVAR( l )
                  WK( j ) = WK( j ) + scalee * WK( n + i )
                  l = l + 1
  110          CONTINUE
            ELSE

!  The IEL-th element has no internal representation.

               DO 120 i = 1, nvarel
                  j = IELVAR( l )
                  WK( j ) = WK( j ) + scalee * GUVALS( k )
                  k = k + 1
                  l = l + 1
  120          CONTINUE
            END IF
  130    CONTINUE

!  Include the contribution from the linear element.

         DO 140 k = ISTADA( ig ), ISTADA( ig1 ) - 1
            j = ICNA( k )
            WK( j ) = WK( j ) + A( k )
  140    CONTINUE

!  The gradient is complete. form the J-th column of the rank-one matrix

         DO 190 l = listvs, listve
            jj = ISVGRP( l )
            j = IWK( jj )
            IF ( j == 0 ) GO TO 190

!  Find the entry in row i of this column.

            DO 180 k = listvs, listve
               ii = ISVGRP( k )
               i = IWK( ii )
               IF ( i == 0 .OR. i > j ) GO TO 180

!  Skip all elements which lie outside a band of width NSEMIB.

               IF ( BAND ) maxsbw = MAX( maxsbw, j - i )
               IF ( j - i > nsemib ) GO TO 180
               hesnew = WK( ii ) * WK( jj ) * g2dash
               IF ( iprint >= 100 ) WRITE( iout, 2090 ) i, j, hesnew
               IF ( NOZERO .AND. hesnew == zero ) GO TO 180

!  Obtain the appropriate storage location in H for the new entry.


!  Case 1: band matrix storage scheme.

               IF ( BAND ) THEN

!  The entry belongs on the diagonal.

                  IF ( i == j ) THEN
                     H( i ) = H( i ) + hesnew

!  The entry belongs off the diagonal.

                  ELSE
                     H( nfree + j - i + nsemib * ( i - 1 ) ) =  &
               H( nfree + j - i + nsemib * ( i - 1 ) ) + hesnew
                  END IF   

!  Case 2: coordinate storage scheme.

               ELSE
                  istart = j
  150             CONTINUE
                  inext = NXTROW( 1, istart )
                  IF ( inext == - 1 ) THEN
                     IF ( newpt > lnxtrw ) GO TO 620

!  The (I,J)-th location is empty. Place the new entry in this location
!  and add another link to the list.

                     nnzh = nnzh + 1
                     IF ( nnzh > lh .OR. nnzh > lih ) GO TO 610
                     IRNH( nnzh ) = i
                     JCNH( nnzh ) = j
                     H( nnzh ) = hesnew
                     NXTROW( 1, istart ) = newpt
                     NXTROW( 2, istart ) = nnzh
                     NXTROW( 1, newpt ) = - 1
                     newpt = newpt + 1
                  ELSE

!  Continue searching the linked list for an entry in row i, column J.

                     IF ( IRNH( NXTROW( 2, istart ) ) == i ) THEN
                        ip = NXTROW( 2, istart )
                        H( ip ) = H( ip ) + hesnew
                     ELSE
                        istart = inext
                        GO TO 150
                     END IF
                  END IF
               END IF
  180       CONTINUE
  190    CONTINUE
  200 CONTINUE

!  Reset the workspace array to zero.

      CALL SETVL( maxsel, WK, 1, zero )

! --------------------------------------------------------
!  Add on the low rank first order terms for the I-th group.
! --------------------------------------------------------

      DO 300 ig = 1, ng
         IF ( iprint >= 100 ) WRITE( iout, 2100 ) ig
         IF ( .NOT. FIXSTR .AND. GSCALE( ig ) == zero ) GO TO 300
         IF ( GXEQX( ig ) ) THEN
            gdash = GSCALE( ig )
         ELSE
            gdash = GSCALE( ig ) * GVALS2( ig )
            IF ( NOZERO .AND. gdash == zero ) GO TO 300
         END IF
         ig1 = ig + 1

!  See if the group has any nonlinear elements.

         DO 290 iell = ISTADG( ig ), ISTADG( ig1 ) - 1
            iel = IELING( iell )
            listvs = ISTAEV( iel )
            listve = ISTAEV( iel + 1 ) - 1
            nvarel = listve - listvs + 1
            ielh = ISTADH( iel )
            ihnext = ielh
            scalee = ESCALE( iell )
            DO 250 l = listvs, listve
               j = IWK( IELVAR( l ) )
               IF ( j /= 0 ) THEN

!  The IEL-th element has an internal representation.
!  Compute the J-th column of the element Hessian matrix.

                  IF ( INTREP( iel ) ) THEN

!  Compute the J-th column of the Hessian.

                     WK( l - listvs + 1 ) = one

!  Find the internal variables.

                     nin = INTVAR( iel + 1 ) - INTVAR( iel )
                     CALL RANGE ( iel, .FALSE., WK( 1 ), &
                              WK( maxsel + 1 ), nvarel, nin, &
                              ITYPEE( iel ), nvarel, nin )

!  Multiply the internal variables by the element Hessian.

                     nn = maxsel + nin
      CALL SETVL( nin, WK( nn + 1 ), 1, zero )

!  Only the upper triangle of the element Hessian is stored.

                     jcolst = ielh - 1
                     DO 220 jcol = 1, nin
                        ijhess = jcolst
                        jcolst = jcolst + jcol
                        wki = WK( maxsel + jcol ) * gdash
                        DO 210 irow = 1, nin
                           IF ( irow <= jcol ) THEN
                              ijhess = ijhess + 1
                           ELSE
                              ijhess = ijhess + irow - 1
                           END IF
                           WK( nn + irow ) = WK( nn + irow ) + &
                                             wki * HUVALS( ijhess )
  210                   CONTINUE
  220                CONTINUE

!  Scatter the product back onto the elemental variables.

                     nnn = nn + nin
                     CALL RANGE ( iel, .TRUE., WK( nn + 1 ), &
                                  WK( nnn + 1 ), nvarel, nin, &
                                  ITYPEE( iel ), nin, nvarel )
                     WK( l - listvs + 1 ) = zero
                  END IF

!  Find the entry in row i of this column.

                  DO 240 k = listvs, l
                     i = IWK( IELVAR( k ) )

!  Skip all elements which lie outside a band of width NSEMIB.

                     IF ( BAND .AND. i /= 0 ) maxsbw =  &
                          MAX( maxsbw, ABS( j - i ) )
                     IF ( ABS( i - j ) <= nsemib .AND. i /= 0 ) THEN

!  Only the upper triangle of the matrix is stored.

                        IF ( i <= j ) THEN
                           ii = i
                           jj = j
                        ELSE
                           ii = j
                           jj = i
                        END IF

!  Obtain the appropriate storage location in H for the new entry.

                        IF ( INTREP( iel ) ) THEN
                           hesnew = scalee * WK( nnn + k - listvs + 1 )
                        ELSE
                           hesnew = scalee * HUVALS( ihnext ) * gdash
                        END IF
                        IF ( iprint >= 100 ) &
                           WRITE( 6, 2080 ) ii, jj, iel, hesnew

!  Case 1: band matrix storage scheme.

                        IF ( BAND ) THEN

!  The entry belongs on the diagonal.

                           IF ( ii == jj ) THEN
                              H( ii ) = H( ii ) + hesnew
                              IF ( k /= l ) H( ii ) = H( ii ) + hesnew

!  The entry belongs off the diagonal.

                           ELSE
                              H( nfree + jj - ii + nsemib * ( ii - 1 ) ) &
 =  H( nfree + jj - ii + nsemib * ( ii - 1 ) ) &
 + hesnew
                           END IF   

!  Case 2: coordinate storage scheme.

                        ELSE
                           IF ( .NOT. NOZERO .OR. hesnew /= ZERO)THEN
                              istart = jj
  230                         CONTINUE
                              inext = NXTROW( 1, istart )
                              IF ( inext == - 1 ) THEN
                                 IF ( newpt > lnxtrw ) GO TO 620

!  The (I,J)-th location is empty. Place the new entry in this location
!  and add another link to the list.

                                 nnzh = nnzh + 1
                                 IF ( nnzh > lh .OR.  &
                                      nnzh > lih ) GO TO 610
                                 IRNH( nnzh ) = ii
                                 JCNH( nnzh ) = jj
                                 H( nnzh ) = hesnew
                                 IF ( k /= l .AND. ii == jj )   &
                                    H( nnzh ) = hesnew + hesnew
                                 NXTROW( 1, istart ) = newpt
                                 NXTROW( 2, istart ) = nnzh
                                 NXTROW( 1, newpt ) = - 1
                                 newpt = newpt + 1
                              ELSE

!  Continue searching the linked list for an entry in row i, column J.

                              IF ( IRNH( NXTROW( 2, istart ) ) &
 == ii ) THEN
                                    ip = NXTROW( 2, istart )
                                    H( ip ) = H( ip ) + hesnew
                                    IF ( k /= l .AND. ii == jj )   &
                                       H( ip ) = H( ip ) + hesnew
                                 ELSE
                                    istart = inext
                                    GO TO 230
                                 END IF
                              END IF
                           END IF
                        END IF
                     END IF
                     ihnext = ihnext + 1
  240             CONTINUE
               ELSE
                  ihnext = ihnext + l - listvs + 1
               END IF
  250       CONTINUE
  290    CONTINUE
  300 CONTINUE

! ----------------------------------------
!  For debugging, print the nonzero values.
! ----------------------------------------

      IF ( iprint >= 10 ) THEN
         IF ( .NOT. BAND ) WRITE( iout, 2000 ) &
 ( IRNH( i ), JCNH( i ), H( i ), i = 1, nnzh )
         IF ( nfree <= 36 ) THEN
            DO 420 i = 1, nfree
               DO 410 j = 1, nfree
                  MATRIX( i, j ) = '  '
  410          CONTINUE
  420       CONTINUE
            IF ( BAND ) THEN
               DO 440 i = 1, nfree
                  IF ( H( i ) /= zero ) MATRIX( i, i ) = ' *'
                  DO 430 j = 1, MIN( nsemib, nfree - i )
                     k = nfree + j + nsemib * ( i - 1 )
                     IF ( H( k ) /= zero ) THEN
                        MATRIX( i + j, i ) = ' *'
                        MATRIX( i, i + j ) = ' *'
                     END IF
  430             CONTINUE
  440          CONTINUE
            ELSE
               DO 450 i = 1, nnzh
                  IF ( IRNH( i ) > nfree ) THEN
                     WRITE( iout, * ) ' ENTRY OUT OF BOUNDS IN ASEMBL ', &
                                      ' ROW NUMBER = ', IRNH( i )
                     STOP
                  END IF
                  IF ( JCNH( i ) > nfree ) THEN
                     WRITE( iout, * ) ' ENTRY OUT OF BOUNDS IN ASEMBL ', &
                                      ' COL NUMBER = ', JCNH( i )
                     STOP
                  END IF
                  MATRIX( IRNH( i ), JCNH( i ) ) = ' *'
                  MATRIX( JCNH( i ), IRNH( i ) ) = ' *'
  450          CONTINUE
            END IF
            WRITE( iout, 2040 ) ( i, i = 1, nfree )
            DO 460 i = 1, nfree
               WRITE( iout, 2050 ) i, ( MATRIX( i, j ), j = 1, nfree )
  460       CONTINUE
         END IF
      END IF
      inform = 0
      RETURN

!  Unsuccessful returns.

  610 CONTINUE
      inform = 1
      IF ( lh <= lih ) THEN
         WRITE( iout, 2010 ) lh
      ELSE
         WRITE( iout, 2030 ) lih
      END IF
      RETURN
  620 CONTINUE
      inform = 2
      WRITE( iout, 2020 ) lnxtrw
      RETURN

!  Non-executable statements

 2000 FORMAT( '    Row  Column    Value        Row  Column    Value ', / &
              ' --- ------ ----- --- ------ ----- ', / &
 ( 2I6, 1P, D24.16, 2I6, 1P, D24.16 ) )
 2010 FORMAT( ' ** Array dimension lh =', I6, ' too small in ASMBL. ' )
 2020 FORMAT( ' ** Array dimension lnxtrw =',I8, ' too small in ASMBL.')
 2030 FORMAT( ' ** Array dimension lih =', I6, ' too small in ASMBL.' )
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
