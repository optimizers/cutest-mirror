! ( Last modified on 23 Dec 2000 at 22:01:38 )

!  ** FOR THE CRAY 2, LINES STARTING 'CDIR$ IVDEP' TELL THE COMPILER TO
!     IGNORE POTENTIAL VECTOR DEPENDENCIES AS THEY ARE KNOWN TO BE O.K.

      SUBROUTINE DHSPRD( n, nn, ng, ntotel, nfree, nvar1, nvar2,  &
                         nbprod, ALLLIN, IVAR, ISTAEV, lstaev, ISTADH,  &
                         lstadh, INTVAR, lntvar, IELING, leling, IELVAR,  &
                         lelvar, ISTAJC, lnstjc, ISELTS, lnelts, ISPTRS,  &
                         lnptrs, IGCOLJ, lngclj, ISLGRP, lnlgrp, ISWKSP,  &
                         lnwksp, ISVGRP, lnvgrp, ISTAGV, lnstgv, IVALJR,  &
                         lnvljr, ITYPEE, litype, nnonnz, INONNZ, lnnnon, &
                         IUSED, lniuse, INONZ2, lnnno2, ISYMMH, maxszh, &
                         P, Q, GVALS2, GVALS3, GRJAC, lgrjac, &
                         GSCALE, ESCALE, lescal, HUVALS, lhuval, &
                         WK, lnwk, WKB, lnwkb, WKC, lnwkc,  &
                         gxeqx, lgxeqx, intrep, lintre, DENSEP, RANGE )
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  EVALUATE Q, THE PRODUCT OF THE HESSIAN OF A GROUP  PARTIALLY
!  SEPARABLE FUNCTION WITH THE VECTOR P.

!  THE NONZERO COMPONENTS OF P HAVE INDICES IVAR( i ), i = nvar1, NVAR2.
!  THE NONZERO COMPONENTS OF THE PRODUCT HAVE INDICES INNONZ( i ),
!  i = 1, NNONNZ.

!  THE ELEMENTS OF THE ARRAY IUSED MUST BE SET TO zero ON ENTRY.
!  THE WORKSPACE ARRAY WK MUST HAVE LENGTH AT LEAST MAX( ng, n + 2 *
!  MAXIMUM NUMBER OF INTERNAL VARIABLES )

!  NICK GOULD, 10TH MAY 1989.
!  FOR CGT PRODUCTIONS.

      INTEGER :: n, nn, ng, ntotel, nfree, nvar1, nvar2, nbprod
      INTEGER :: nnonnz, maxszh, lgxeqx, lintre, litype
      INTEGER :: lstaev, lstadh, lntvar, leling, lelvar
      INTEGER :: lnstjc, lnelts, lnptrs, lngclj, lnlgrp, lnwksp
      INTEGER :: lnvgrp, lnstgv, lnvljr, lnnnon, lniuse, lnnno2
      INTEGER :: lgrjac, lescal, lhuval, lnwk, lnwkb, lnwkc
      LOGICAL :: ALLLIN, DENSEP
      INTEGER :: IVAR( n ), ISTAEV( lstaev ), ISTADH( lstadh )
      INTEGER :: ISTAJC( lnstjc ), ISELTS( lnelts )
      INTEGER :: IGCOLJ( lngclj ), ISLGRP( lnlgrp )
      INTEGER :: ISWKSP( lnwksp ), ISPTRS( lnptrs )
      INTEGER :: INTVAR( lntvar ), IELING( leling )
      INTEGER :: IELVAR( lelvar ), ISVGRP( lnvgrp )
      INTEGER :: ISTAGV( lnstgv ), IVALJR( lnvljr )
      INTEGER :: INONNZ( lnnnon ), IUSED( lniuse )
      INTEGER :: ITYPEE( litype )
      INTEGER :: INONZ2( lnnno2 ), ISYMMH( maxszh, maxszh )
      REAL ( KIND = wp ) :: P( n ), GVALS2( ng ), GVALS3( ng ), &
                       GSCALE( ng ), ESCALE( lescal ), &
                       WK( lnwk ), WKB( lnwkb ), WKC( lnwkc ), &
                       Q( n ), GRJAC( lgrjac ), HUVALS( lhuval )
      LOGICAL :: GXEQX( lgxeqx ), INTREP( lintre )
      EXTERNAL :: RANGE 

!  LOCAL VARIABLES.

      INTEGER :: i, iel, ig, ipt, j, irow, jcol, ijhess, lthvar
      INTEGER :: iell, ii, k, l, ll, nin, nvarel, ielhst, nnonz2
      REAL ( KIND = wp ) :: zero, pi, gi
      LOGICAL :: NULLWK
!D    EXTERNAL         SETVL
      INTRINSIC        ABS

!  COMMON VARIABLES.

      REAL ( KIND = wp ) :: epsmch, epsneg, tiny, big
!D    COMMON / DMACHN / epsmch, epsneg, tiny, big
!D    SAVE   / DMACHN /

!  SET CONSTANT REAL PARAMETERS.

      PARAMETER ( zero = 0.0_wp )

! ======================= RANK-ONE TERMS ==========================

!  IF THE IG-TH GROUP is NON-TRIVIAL, FORM THE PRODUCT OF P WITH THE
!  SUM OF RANK-ONE FIRST ORDER TERMS, A(TRANS) * GVALS3 * A. A is
!  STORED BY BOTH ROWS AND COLUMNS. FOR MAXIMUM EFFICIENCY,
!  THE PRODUCT is FORMED IN DIFFERENT WAYS IF P is SPARSE OR DENSE.

! -----------------  CASE 1. P is NOT SPARSE -----------------------

      IF ( DENSEP ) THEN

!  INITIALIZE WK AND Q AS ZERO.

      CALL SETVL( ng, WK, 1, zero )
      CALL SETVL( n,  Q,  1, zero )

!  FORM THE MATRIX-VECTOR PRODUCT WK = A * P, USING THE COLUMN-WISE
!  STORAGE OF A.

         DO 20 j = nvar1, nvar2
            i = IVAR( j )
            pi = P( i )
!DIR$ IVDEP
            DO 10 k = ISTAJC( i ), ISTAJC( i + 1 ) - 1
               l = IGCOLJ( k )
               WK( l ) = WK( l ) + pi * GRJAC( k )
   10       CONTINUE
   20    CONTINUE

!  MULTIPLY WK BY THE DIAGONAL MATRIX GVALS3.

         DO 30 ig = 1, ng
            IF ( GXEQX( ig ) ) THEN
               WK( ig ) = GSCALE( ig ) * WK( ig )
            ELSE
               WK( ig ) = GSCALE( ig ) * WK( ig ) * GVALS3( ig )
            END IF
   30    CONTINUE

!  FORM THE MATRIX-VECTOR PRODUCT Q = A(TRANS) * WK, ONCE
!  AGAIN USING THE COLUMN-WISE STORAGE OF A.

         nnonnz = 0
         DO 50 j = 1, nfree
            i = IVAR( j )
            pi = zero
!DIR$ IVDEP
            DO  40 k = ISTAJC( i ), ISTAJC( i + 1 ) - 1
               pi = pi + WK( IGCOLJ( k ) ) * GRJAC( k )
   40       CONTINUE
            Q( i ) = pi
   50    CONTINUE
      ELSE

! ------------------- CASE 2. P is SPARSE --------------------------

         nnonz2 = 0

!  FORM THE MATRIX-VECTOR PRODUCT WK = A * P, USING THE COLUMN-WISE
!  STORAGE OF A. KEEP TRACK OF THE NONZERO COMPONENTS OF WK IN INONZ2.
!  ONLY STORE COMPONENTS CORRESPONDING TO NON TRIVIAL GROUP .

         DO 120 j = nvar1, nvar2
            i = IVAR( j )
            pi = P( i )
!DIR$ IVDEP
            DO 110 k = ISTAJC( i ), ISTAJC( i + 1 ) - 1
               ig = IGCOLJ( k )
               IF ( IUSED( ig ) == 0 ) THEN
                  WK( ig ) = pi * GRJAC( k )
                  IUSED( ig ) = 1
                  nnonz2 = nnonz2 + 1
                  INONZ2( nnonz2 ) = ig
               ELSE
                  WK( ig ) = WK( ig ) + pi * GRJAC( k )
               END IF
  110       CONTINUE
  120    CONTINUE

!  RESET IUSED TO ZERO.

         DO 130 i = 1, nnonz2
            IUSED( INONZ2( i ) ) = 0
  130    CONTINUE

!  FORM THE MATRIX-VECTOR PRODUCT Q = A(TRANS) * WK, USING THE
!  ROW-WISE STORAGE OF A.

         nnonnz = 0
         DO 160 j = 1, nnonz2
            ig = INONZ2( j )
            IF ( .NOT. GXEQX( ig ) ) THEN

!  IF GROUP ig is NON TRIVIAL, THERE ARE CONTRIBUTIONS FROM ITS
!  RANK-ONE TERM.

               pi = GSCALE( ig ) * GVALS3( ig ) * WK( ig )
!DIR$ IVDEP
               DO 150 k = ISTAGV( ig ), ISTAGV( ig + 1 ) - 1
                  l = ISVGRP( k )

!  IF Q HAS A NONZERO IN POSITION l, STORE ITS INDEX IN INONNZ.

                  IF ( IUSED( l ) == 0 ) THEN
                     Q( l ) = pi * GRJAC( IVALJR( k ) )
                     IUSED( l ) = 1
                     nnonnz = nnonnz + 1
                     INONNZ( nnonnz ) = l
                  ELSE
                     Q( l ) = Q( l ) + pi * GRJAC( IVALJR( k ) )
                  END IF
  150          CONTINUE
            END IF
  160    CONTINUE
      END IF
      IF ( .NOT. ALLLIN ) THEN

! ======================= SECOND-ORDER TERMS =======================

!  NOW CONSIDER THE PRODUCT OF P WITH THE SECOND ORDER TERMS (THAT IS.
! (2ND DERIVATIVES OF THE ELEMENTS). AGAIN, FOR MAXIMUM EFFICIENCY,
!  THE PRODUCT is FORMED IN DIFFERENT WAYS IF P is SPARSE OR DENSE.

! --------------------- CASE 1. P is NOT SPARSE ---------------------

         IF ( DENSEP ) THEN
            DO 280 iell = 1, ntotel
               iel = IELING( iell )
               ig = ISLGRP( iell )
               nvarel = ISTAEV( iel + 1 ) - ISTAEV( iel )
               IF ( GXEQX( ig ) ) THEN
                  gi = GSCALE( ig ) * ESCALE( iell )
               ELSE
                  gi = GSCALE( ig ) * ESCALE( iell ) * GVALS2( ig ) 
               END IF
               ISWKSP( iel ) = nbprod
               IF ( INTREP( iel ) ) THEN

!  THE IEL-TH ELEMENT HESSIAN HAS AN INTERNAL REPRESENTATION.
!  COPY THE ELEMENTAL VARIABLES INTO WK.

                  NULLWK = .TRUE.
                  ll = ISTAEV( iel )
!DIR$ IVDEP
                  DO 210 ii = 1, nvarel
                     pi = P( IELVAR( ll ) )
                     WK( ii ) = pi
                     IF ( pi /= zero ) NULLWK = .FALSE.
                     ll = ll + 1
  210             CONTINUE

!  SKIP THE ELEMENT IF WK is NULL.

                  IF ( NULLWK ) GO TO 280

!  FIND THE INTERNAL VARIABLES, WKB.

                  nin = INTVAR( iel + 1 ) - INTVAR( iel )
                  CALL RANGE ( iel, .FALSE., WK, WKB, nvarel, nin,  &
                               ITYPEE( iel ), nvarel, nin )

!  MULTIPLY THE INTERNAL VARIABLES BY THE ELEMENT HESSIAN AND PUT THE
!  PRODUCT IN WKC. CONSIDER THE FIRST COLUMN OF THE ELEMENT HESSIAN.

                  ielhst = ISTADH( iel )
                  pi = gi * WKB( 1 )
!DIR$ IVDEP
                  DO 220 irow = 1, nin
                     ijhess = ISYMMH( 1, irow ) + ielhst
                     WKC( irow ) = pi * HUVALS( ijhess )
  220             CONTINUE

!  NOW CONSIDER THE REMAINING COLUMNS OF THE ELEMENT HESSIAN.

                  DO 240 jcol = 2, nin
                     pi = gi * WKB( jcol )
                     IF ( pi /= zero ) THEN
!DIR$ IVDEP
                        DO 230 irow = 1, nin
                           ijhess = ISYMMH( jcol, irow ) + ielhst
                           WKC( irow ) = WKC( irow ) + pi * &
                                         HUVALS( ijhess )
  230                   CONTINUE
                     END IF
  240             CONTINUE

!  SCATTER THE PRODUCT BACK ONTO THE ELEMENTAL VARIABLES, WK.

                  CALL RANGE ( iel, .TRUE., WKC, WK, nvarel, nin,  &
                               ITYPEE( iel ), nin, nvarel )

!  ADD THE SCATTERED PRODUCT TO Q.

                  ll = ISTAEV( iel )
!DIR$ IVDEP
                  DO 250 ii = 1, nvarel
                     l = IELVAR( ll )
                     Q( l ) = Q( l ) + WK( ii )
                     ll = ll + 1
  250             CONTINUE
               ELSE

!  THE IEL-TH ELEMENT HESSIAN HAS NO INTERNAL REPRESENTATION.

                  lthvar = ISTAEV( iel ) - 1
                  ielhst = ISTADH( iel )
                  DO 270 jcol = 1, nvarel
                     pi = gi * P( IELVAR( lthvar + jcol ) )
                     IF ( pi /= zero ) THEN
!DIR$ IVDEP
                        DO 260 irow = 1, nvarel
                           ijhess = ISYMMH( jcol, irow ) + ielhst
                           l = IELVAR( lthvar + irow )
                           Q( l ) = Q( l ) + pi * HUVALS( ijhess )
  260                   CONTINUE
                     END IF
  270             CONTINUE
               END IF
  280       CONTINUE
         ELSE

! -------------------- CASE 2. P is SPARSE ------------------------

            DO 380 j = nvar1, nvar2

!  CONSIDER EACH NONZERO COMPONENT OF P SEPARATELY.

               i = IVAR( j )
               ipt = ISPTRS( i )
               IF ( ipt >= 0 ) THEN

!  THE INDEX OF THE I-TH COMPONENT LIES IN THE IEL-TH NONLINEAR ELEMENT.

                  iell = ISELTS( i )
  300             CONTINUE

!  CHECK TO ENSURE THAT THE IEL-TH ELEMENT HAS NOT ALREADY BEEN USED.

                  IF ( ISWKSP( iell ) < nbprod ) THEN
                     ISWKSP( iell ) = nbprod
                     iel = IELING( iell )
                     nvarel = ISTAEV( iel + 1 ) - ISTAEV( iel )
                     ig = ISLGRP( iell )
                     IF ( GXEQX( ig ) ) THEN
                        gi = GSCALE( ig ) * ESCALE( iell )
                     ELSE
                        gi = GSCALE( ig ) * ESCALE( iell ) * GVALS2( IG) 
                     END IF
                     IF ( INTREP( iel ) ) THEN

!  THE IEL-TH ELEMENT HESSIAN HAS AN INTERNAL REPRESENTATION.
!  COPY THE ELEMENTAL VARIABLES INTO WK.

                        ll = ISTAEV( iel )
!DIR$ IVDEP
                        DO 310 ii = 1, nvarel
                           WK( ii ) = P( IELVAR( ll ) )
                           ll = ll + 1
  310                   CONTINUE

!  FIND THE INTERNAL VARIABLES.

                        nin = INTVAR( iel + 1 ) - INTVAR( iel )
                        CALL RANGE ( iel, .FALSE., WK, WKB, nvarel, nin, &
                                    ITYPEE( iel ), nvarel, nin )

!  MULTIPLY THE INTERNAL VARIABLES BY THE ELEMENT HESSIAN AND PUT THE 
!  PRODUCT IN WKB. CONSIDER THE FIRST COLUMN OF THE ELEMENT HESSIAN.

                        ielhst = ISTADH( iel )
                        pi = gi * WKB( 1 )
!DIR$ IVDEP
                        DO 320 irow = 1, nin
                           ijhess = ISYMMH( 1, irow ) + ielhst
                           WKC( irow ) = pi * HUVALS( ijhess )
  320                   CONTINUE

!  NOW CONSIDER THE REMAINING COLUMNS OF THE ELEMENT HESSIAN.

                        DO 340 jcol = 2, nin
                           pi = gi * WKB( jcol )
                           IF ( pi /= zero ) THEN
!DIR$ IVDEP
                              DO 330 irow = 1, nin
                                 ijhess = ISYMMH( jcol, irow ) + ielhst
                                 WKC( irow ) = WKC( irow ) + pi * &
                                             HUVALS( ijhess )
  330                         CONTINUE
                           END IF
  340                   CONTINUE

!  SCATTER THE PRODUCT BACK ONTO THE ELEMENTAL VARIABLES, WK.

                        CALL RANGE ( iel, .TRUE., WKC, WK, nvarel, nin, &
                                     ITYPEE( iel ), nin, nvarel )

!  ADD THE SCATTERED PRODUCT TO Q.

                        ll = ISTAEV( iel )
!DIR$ IVDEP
                        DO 350 ii = 1, nvarel
                           l = IELVAR( ll )

!  IF Q HAS A NONZERO IN POSITION l, STORE ITS INDEX IN INONNZ.

                           IF ( ABS( WK( ii ) ) > tiny ) THEN
                              IF ( IUSED( l ) == 0 ) THEN
                                 Q( l ) = WK( ii )
                                 IUSED( l ) = 1
                                 nnonnz = nnonnz + 1
                                 INONNZ( nnonnz ) = l
                              ELSE
                                 Q( l ) = Q( l ) + WK( ii )
                              END IF
                           END IF
                           ll = ll + 1
  350                   CONTINUE
                     ELSE

!  THE IEL-TH ELEMENT HESSIAN HAS NO INTERNAL REPRESENTATION.

                        lthvar = ISTAEV( iel ) - 1
                        ielhst = ISTADH( iel )
                        DO 370 jcol = 1, nvarel
                           pi = gi * P( IELVAR( lthvar + jcol ) )
                           IF ( pi /= zero ) THEN
!DIR$ IVDEP
                              DO 360 irow = 1, nvarel
                                 ijhess = ISYMMH( jcol, irow ) + ielhst

!  IF Q HAS A NONZERO IN POSITION l, STORE ITS INDEX IN INONNZ.

                                 IF ( ABS( HUVALS( ijhess ) ) &
                                    > tiny ) THEN
                                    l = IELVAR( lthvar + irow )
                                    IF ( IUSED( l ) == 0 ) THEN
                                       Q( l ) = pi * HUVALS( ijhess )
                                       IUSED( l ) = 1
                                       nnonnz = nnonnz + 1
                                       INONNZ( nnonnz ) = l
                                    ELSE
                                       Q( l ) = Q( l ) + pi * &
                                                HUVALS( ijhess )
                                    END IF
                                 END IF
  360                         CONTINUE
                           END IF
  370                   CONTINUE
                     END IF
                  END IF

!  CHECK TO SEE IF THERE ARE ANY FURTHER ELEMENTS WHOSE VARIABLES
!  INCLUDE THE I-TH VARIABLE.

                  IF ( ipt > 0 ) THEN
                     iell = ISELTS( ipt )
                     ipt = ISPTRS( ipt )
                     GO TO 300
                  END IF
               END IF
  380       CONTINUE
         END IF
      END IF

! ==================== THE PRODUCT is COMPLETE =======================

!  RESET IUSED TO ZERO.

      IF ( .NOT. DENSEP ) THEN
         DO 390 i = 1, nnonnz
            IUSED( INONNZ( i ) ) = 0
  390    CONTINUE
      END IF
      RETURN

!  END OF HSPRD.

      END
