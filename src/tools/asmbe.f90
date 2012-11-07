! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE ASMBE( n, ng, maxsel, ISTADH, lstadh, ICNA, licna, &
                         ISTADA, lstada, INTVAR, lntvar, IELVAR, lelvar, &
                         IELING, leling, ISTADG, lstadg, ISTAEV, lstaev, &
                         ISTAGV, lnstgv, ISVGRP, lnvgrp, IWK, liwk,  &
                         A, la, GUVALS, lnguvl, HUVALS, lnhuvl, GVALS2,  &
                         GVALS3, GSCALE, ESCALE, lescal, WK, lwk,  &
                         gxeqx, lgxeqx, intrep, lintre, ITYPEE, litype, &
                         RANGE, ne, IRNHI, lirnhi, IPRNHI, HI,  &
                         lhi, IPRHI, BYROWS, iprint, iout, INFORM)
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  ********************************************************************

!  Assemble the second derivative matrix of a groups partially
!  separable function into finite-element format

!  Nick Gould, November 25th 1994, for CGT Productions.

!  ********************************************************************

      INTEGER :: n, ng, maxsel, lirnhi, lhi, ne
      INTEGER :: la, litype, iout, inform, iprint
      INTEGER :: lstadh, licna, lstada, lntvar, lelvar, leling
      INTEGER :: lstadg, lstaev, lnstgv, lnvgrp, liwk
      INTEGER :: lnguvl, lnhuvl, lescal, lwk, lgxeqx, lintre
      LOGICAL :: BYROWS
      INTEGER :: ISTADH( lstadh ), ICNA ( licna )
      INTEGER :: ISTADA( lstada ), INTVAR( lntvar )
      INTEGER :: IELVAR( lelvar ), IELING( leling )
      INTEGER :: ISTADG( lstadg ), ISTAEV( lstaev )
      INTEGER :: ISTAGV( lnstgv ), ISVGRP( lnvgrp )
      INTEGER :: IRNHI ( lirnhi ), IWK ( liwk )
      INTEGER :: IPRNHI( ng + 1 ), IPRHI ( ng + 1 )
      INTEGER :: ITYPEE( litype )
      LOGICAL :: GXEQX( lgxeqx ),  INTREP( lintre )
      REAL ( KIND = wp ) :: A(     la ),      GVALS2( ng ),     GVALS3( ng ), &
                       GUVALS( lnguvl ), HUVALS( lnhuvl ), GSCALE( ng ), &
                       ESCALE( lescal ), HI ( lhi ), WK(     lwk )
      EXTERNAL :: RANGE 

!  Local variables.

      INTEGER :: i, ii, j, ihi, jj, k, l, nn, ig, nnn
      INTEGER :: nin,    iell,   iel, ihnext, nvarg, nsizee
      INTEGER :: nvarel, ig1,    listvs, listve, ielh
      INTEGER :: ijhess, irow,   jcol,   jcolst
      REAL ( KIND = wp ) :: one,    zero,   wki,    hesnew, gdash,  g2dash, &
                       scalee
!D    EXTERNAL         SETVL, SETVI

!  Set constant real parameters.

      PARAMETER ( zero = 0.0_wp, one = 1.0_wp )

! ------------------------------------------------------
!  Form the rank-one second order term for the i-th group.
! ------------------------------------------------------

      ne = 0
      IPRNHI( 1 ) = 1
      IPRHI ( 1 ) = 1
      DO 200 ig = 1, ng
         IF ( iprint >= 100 ) WRITE( iout, 2070 ) ig
         ig1 = ig + 1
         IF ( GXEQX( ig ) ) THEN
            g2dash = zero
         ELSE
            g2dash = GSCALE( ig ) * GVALS3( ig )
            IF ( iprint >= 100 )WRITE(6,*) ' GVALS3(IG) ', GVALS3(IG)
         END IF

!  Ignore linear groups

         IF ( ISTADG( ig ) >= ISTADG( ig1 ) .AND. GXEQX( ig ) )  &
            GO TO 200

!  The group is nonlinear

         ne = ne + 1
         listvs = ISTAGV( ig )
         listve = ISTAGV( ig1 ) - 1
         nvarg = listve - listvs + 1

!  Set the starting addresses for the integer and real arrays for
!  group ig + 1

         IPRNHI( ne + 1 ) = IPRNHI( ne ) + nvarg
         IF ( IPRNHI( ne + 1 ) > lirnhi ) THEN
            WRITE( iout, 2030 ) lirnhi
            GO TO 610
         END IF
         nsizee = ( nvarg * ( nvarg + 1 ) ) / 2
         IPRHI( ne + 1 ) = IPRHI( ne ) + nsizee
         IF ( IPRHI( ne + 1 ) >= lhi ) THEN
            WRITE( iout, 2040 ) lhi, IPRHI( ne + 1 )
            GO TO 610
         END IF

!  Record the row indices involved in super-element ne

         IF ( GXEQX( ig ) .OR. g2dash == zero ) THEN 
      CALL SETVL( nsizee, HI( IPRHI( ne ) ), 1, zero )
            GO TO 200
         END IF
         k = IPRNHI( ne )
         DO 110 l = listvs, listve
            IRNHI( k ) = ISVGRP( l )
            k = k + 1
  110    CONTINUE   

!  Form the gradient of the ig-th group.

      CALL SETVI( nvarg, WK, ISVGRP( listvs ), zero )

!  Consider any nonlinear elements for the group.

         DO 160 iell = ISTADG( ig ), ISTADG( ig1 ) - 1
            iel = IELING( iell )
            k = INTVAR( iel )
            l = ISTAEV( iel )
            nvarel = ISTAEV( iel + 1 ) - l
            scalee = ESCALE( iell )
            IF ( INTREP( iel ) ) THEN

!  The IEL-th element has an internal representation.

               nin = INTVAR( iel + 1 ) - k
               CALL RANGE ( iel, .TRUE., GUVALS( k ), &
                            WK( n + 1 ), nvarel, nin, &
                            ITYPEE( iel ), nin, nvarel )
               DO 140 i = 1, nvarel
                  j = IELVAR( l )
                  WK( j ) = WK( j ) + scalee * WK( n + i )
                  l = l + 1
  140          CONTINUE
            ELSE

!  The IEL-th element has no internal representation.

               DO 150 i = 1, nvarel
                  j = IELVAR( l )
                  WK( j ) = WK( j ) + scalee * GUVALS( k )
                  k = k + 1
                  l = l + 1
  150          CONTINUE
            END IF
  160    CONTINUE

!  Include the contribution from the linear element.

         DO 170 k = ISTADA( ig ), ISTADA( ig1 ) - 1
            j = ICNA( k )
            WK( j ) = WK( j ) + A( k )
  170    CONTINUE

!  The gradient is complete. Form the J-th column of the rank-one matrix

         DO 190 l = listvs, listve
            jj = ISVGRP( l )
            j = l - listvs + 1

!  Find the entry in row i of this column.

            DO 180 k = listvs, l
               ii = ISVGRP( k )
               i = k - listvs + 1
               IF ( BYROWS ) THEN
                  ihi = IPRHI( ne ) - 1 + nvarg * ( i - 1 ) -  &
 ( ( i - 1 ) * i ) / 2 + j
               ELSE
                  ihi = IPRHI( ne ) - 1 + i + ( j * ( j - 1 ) ) / 2
               END IF
               HI( ihi ) = WK( ii ) * WK( jj ) * g2dash
               IF ( iprint >= 100 )  &
                  WRITE( iout, 2090 ) ii, jj, HI( ihi )
  180       CONTINUE
  190    CONTINUE
  200 CONTINUE

!  Reset the workspace array to zero.

      CALL SETVL( maxsel, WK, 1, zero )

! --------------------------------------------------------
!  Add on the low rank first order terms for the I-th group.
! --------------------------------------------------------

      ne = 0
      DO 300 ig = 1, ng
         ig1 = ig + 1

!  Once again, ignore linear groups

         IF ( ISTADG( ig ) >= ISTADG( ig1 ) .AND. GXEQX( ig ) )  &
            GO TO 300

!  The group is nonlinear

         ne = ne + 1
         IF ( iprint >= 100 ) WRITE( iout, 2100 ) ig
         IF ( GXEQX( ig ) ) THEN
            gdash = GSCALE( ig )
         ELSE
            gdash = GSCALE( ig ) * GVALS2( ig )
            IF ( iprint >= 100 )WRITE(6,*) ' GVALS2(IG) ', GVALS2(IG)
         END IF
         IF ( gdash == zero ) GO TO 300

!  Map the problem variables to the elemental variables.

         nvarg = IPRNHI( ne + 1 ) - IPRNHI( ne )
         DO 210 i = IPRNHI( ne ), IPRNHI( ne + 1 ) - 1
            IWK( IRNHI( i ) ) = i + 1 - IPRNHI( ne )
  210    CONTINUE   

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
                  DO 230 jcol = 1, nin
                     ijhess = jcolst
                     jcolst = jcolst + jcol
                     wki = WK( maxsel + jcol ) * gdash
                     DO 220 irow = 1, nin
                        IF ( irow <= jcol ) THEN
                           ijhess = ijhess + 1
                        ELSE
                           ijhess = ijhess + irow - 1
                        END IF
                        WK( nn + irow ) = WK( nn + irow ) + &
                                          wki * HUVALS( ijhess )
  220                CONTINUE
  230             CONTINUE

!  Scatter the product back onto the elemental variables.

                  nnn = nn + nin
                  CALL RANGE ( iel, .TRUE., WK( nn + 1 ), &
                               WK( nnn + 1 ), nvarel, nin, &
                               ITYPEE( iel ), nin, nvarel )
                  WK( l - listvs + 1 ) = zero

!  Find the entry in row i of this column.

               END IF
               DO 240 k = listvs, l
                  i = IWK( IELVAR( k ) )

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
                  IF ( BYROWS ) THEN
                     ihi = IPRHI( ne ) - 1 + nvarg * ( ii - 1 ) -  &
 ( ( ii - 1 ) * ii ) / 2 + jj
                  ELSE
                     ihi = IPRHI( ne ) - 1 + ii + &
 ( jj * ( jj - 1 ) ) / 2
                  END IF
                  HI( ihi ) = HI( ihi ) + hesnew
                  IF ( k /= l .AND. ii == jj ) &
                       HI( ihi ) = HI( ihi ) + hesnew
                  ihnext = ihnext + 1
  240          CONTINUE
  250       CONTINUE
  290    CONTINUE
  300 CONTINUE

! ----------------------------------------
!  For debugging, print the nonzero values.
! ----------------------------------------

      IF ( iprint >= 10 ) THEN
         DO 410 ig = 1, ne
            WRITE( iout, 2000 ) ig
            WRITE( iout, 2010 ) &
 ( IRNHI( i ), i = IPRNHI( ig ), IPRNHI( ig + 1 ) - 1 )
            WRITE( iout, 2020 ) &
 ( HI( i ), i = IPRHI( ig ), IPRHI( ig + 1 ) - 1 )
  410    CONTINUE   
      END IF
      inform = 0
      RETURN

!  Unsuccessful returns.

  610 CONTINUE
      inform = 1
      RETURN

!  Non-executable statements

 2000 FORMAT( ' Super-element ', I10 )   
 2010 FORMAT( ' Super-element variables     ', 8I7, /, ( 11I7 ) )
 2020 FORMAT( ' Nonzeros   ', 1P, 6D12.4, /, ( 7D12.4 ) )
 2030 FORMAT( ' ** Array dimension lirnhi = ',I8,' too small in ASMBE.')
 2040 FORMAT( ' ** Array dimension lhi = ', I8,' too small in ASMBE.  &
           Increase to at least ', I8 )
 2070 FORMAT( ' Group ', I5, ' rank-one terms ' )
 2080 FORMAT( ' Row ', I6, ' Column ', I6, ' used from element ', I6, &
              ' value = ', 1P, D24.16 )
 2090 FORMAT( ' Row ', I6, ' column ', I6, ' used. Value = ',1P,D24.16)
 2100 FORMAT( ' Group ', I5, ' second-order terms ' )

!  End of ASMBE.

      END
