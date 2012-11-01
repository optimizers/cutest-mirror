! ( Last modified on 23 Dec 2000 at 22:01:38 )

!  ** FOR THE CRAY 2, LINES STARTING 'CDIR$ IVDEP' TELL THE COMPILER TO
!     IGNORE POTENTIAL VECTOR DEPENDENCIES AS THEY ARE KNOWN TO BE O.K.

      SUBROUTINE DHSPRD( N, NN, NG, NGEL, NFREE, NVAR1, NVAR2,  &
                         NBPROD, ALLLIN, IVAR, ISTAEV, LSTAEV, ISTADH,  &
                         LSTADH, INTVAR, LNTVAR, IELING, LELING, IELVAR,  &
                         LELVAR, ISTAJC, LNSTJC, ISELTS, LNELTS, ISPTRS,  &
                         LNPTRS, IGCOLJ, LNGCLJ, ISLGRP, LNLGRP, ISWKSP,  &
                         LNWKSP, ISVGRP, LNVGRP, ISTAGV, LNSTGV, IVALJR,  &
                         LNVLJR, ITYPEE, LITYPE, NNONNZ, INONNZ, LNNNON, &
                         IUSED, LNIUSE, INONZ2, LNNNO2, ISYMMH, MAXSZH, &
                         P, Q, GVALS2, GVALS3, GRJAC, LGRJAC, &
                         GSCALE, ESCALE, LESCAL, HUVALS, LHUVAL, &
                         WK, LNWK, WKB, LNWKB, WKC, LNWKC,  &
                         GXEQX, LGXEQX, INTREP, LINTRE, DENSEP, RANGE )
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  EVALUATE Q, THE PRODUCT OF THE HESSIAN OF A GROUP  PARTIALLY
!  SEPARABLE FUNCTION WITH THE VECTOR P.

!  THE NONZERO COMPONENTS OF P HAVE INDICES IVAR( I ), I = NVAR1, NVAR2.
!  THE NONZERO COMPONENTS OF THE PRODUCT HAVE INDICES INNONZ( I ),
!  I = 1, NNONNZ.

!  THE ELEMENTS OF THE ARRAY IUSED MUST BE SET TO ZERO ON ENTRY.
!  THE WORKSPACE ARRAY WK MUST HAVE LENGTH AT LEAST MAX( NG, N + 2 *
!  MAXIMUM NUMBER OF INTERNAL VARIABLES )

!  NICK GOULD, 10TH MAY 1989.
!  FOR CGT PRODUCTIONS.

      INTEGER :: N, NN, NG, NGEL, NFREE, NVAR1, NVAR2, NBPROD
      INTEGER :: NNONNZ, MAXSZH, LGXEQX, LINTRE, LITYPE
      INTEGER :: LSTAEV, LSTADH, LNTVAR, LELING, LELVAR
      INTEGER :: LNSTJC, LNELTS, LNPTRS, LNGCLJ, LNLGRP, LNWKSP
      INTEGER :: LNVGRP, LNSTGV, LNVLJR, LNNNON, LNIUSE, LNNNO2
      INTEGER :: LGRJAC, LESCAL, LHUVAL, LNWK, LNWKB, LNWKC
      LOGICAL :: ALLLIN, DENSEP
      INTEGER :: IVAR( N ), ISTAEV( LSTAEV ), ISTADH( LSTADH )
      INTEGER :: ISTAJC( LNSTJC ), ISELTS( LNELTS )
      INTEGER :: IGCOLJ( LNGCLJ ), ISLGRP( LNLGRP )
      INTEGER :: ISWKSP( LNWKSP ), ISPTRS( LNPTRS )
      INTEGER :: INTVAR( LNTVAR ), IELING( LELING )
      INTEGER :: IELVAR( LELVAR ), ISVGRP( LNVGRP )
      INTEGER :: ISTAGV( LNSTGV ), IVALJR( LNVLJR )
      INTEGER :: INONNZ( LNNNON ), IUSED( LNIUSE )
      INTEGER :: ITYPEE( LITYPE )
      INTEGER :: INONZ2( LNNNO2 ), ISYMMH( MAXSZH, MAXSZH )
      REAL ( KIND = wp ) :: P( N ), GVALS2( NG ), GVALS3( NG ), &
                       GSCALE( NG ), ESCALE( LESCAL ), &
                       WK( LNWK ), WKB( LNWKB ), WKC( LNWKC ), &
                       Q( N ), GRJAC( LGRJAC ), HUVALS( LHUVAL )
      LOGICAL :: GXEQX( LGXEQX ), INTREP( LINTRE )
      EXTERNAL :: RANGE 

!  LOCAL VARIABLES.

      INTEGER :: I, IEL, IG, IPT, J, IROW, JCOL, IJHESS, LTHVAR
      INTEGER :: IELL, II, K, L, LL, NIN, NVAREL, IELHST, NNONZ2
      REAL ( KIND = wp ) :: ZERO, PI, GI
      LOGICAL :: NULLWK
!D    EXTERNAL         DSETVL
      INTRINSIC        ABS

!  COMMON VARIABLES.

      REAL ( KIND = wp ) :: EPSMCH, EPSNEG, TINY, BIG
!D    COMMON / DMACHN / EPSMCH, EPSNEG, TINY, BIG
!D    SAVE   / DMACHN /

!  SET CONSTANT REAL PARAMETERS.

      PARAMETER ( ZERO = 0.0_wp )

! ======================= RANK-ONE TERMS ==========================

!  IF THE IG-TH GROUP IS NON-TRIVIAL, FORM THE PRODUCT OF P WITH THE
!  SUM OF RANK-ONE FIRST ORDER TERMS, A(TRANS) * GVALS3 * A. A IS
!  STORED BY BOTH ROWS AND COLUMNS. FOR MAXIMUM EFFICIENCY,
!  THE PRODUCT IS FORMED IN DIFFERENT WAYS IF P IS SPARSE OR DENSE.

! -----------------  CASE 1. P IS NOT SPARSE -----------------------

      IF ( DENSEP ) THEN

!  INITIALIZE WK AND Q AS ZERO.

      CALL DSETVL( NG, WK, 1, ZERO )
      CALL DSETVL( N,  Q,  1, ZERO )

!  FORM THE MATRIX-VECTOR PRODUCT WK = A * P, USING THE COLUMN-WISE
!  STORAGE OF A.

         DO 20 J = NVAR1, NVAR2
            I = IVAR( J )
            PI = P( I )
!DIR$ IVDEP
            DO 10 K = ISTAJC( I ), ISTAJC( I + 1 ) - 1
               L = IGCOLJ( K )
               WK( L ) = WK( L ) + PI * GRJAC( K )
   10       CONTINUE
   20    CONTINUE

!  MULTIPLY WK BY THE DIAGONAL MATRIX GVALS3.

         DO 30 IG = 1, NG
            IF ( GXEQX( IG ) ) THEN
               WK( IG ) = GSCALE( IG ) * WK( IG )
            ELSE
               WK( IG ) = GSCALE( IG ) * WK( IG ) * GVALS3( IG )
            END IF
   30    CONTINUE

!  FORM THE MATRIX-VECTOR PRODUCT Q = A(TRANS) * WK, ONCE
!  AGAIN USING THE COLUMN-WISE STORAGE OF A.

         NNONNZ = 0
         DO 50 J = 1, NFREE
            I = IVAR( J )
            PI = ZERO
!DIR$ IVDEP
            DO  40 K = ISTAJC( I ), ISTAJC( I + 1 ) - 1
               PI = PI + WK( IGCOLJ( K ) ) * GRJAC( K )
   40       CONTINUE
            Q( I ) = PI
   50    CONTINUE
      ELSE

! ------------------- CASE 2. P IS SPARSE --------------------------

         NNONZ2 = 0

!  FORM THE MATRIX-VECTOR PRODUCT WK = A * P, USING THE COLUMN-WISE
!  STORAGE OF A. KEEP TRACK OF THE NONZERO COMPONENTS OF WK IN INONZ2.
!  ONLY STORE COMPONENTS CORRESPONDING TO NON TRIVIAL GROUP .

         DO 120 J = NVAR1, NVAR2
            I = IVAR( J )
            PI = P( I )
!DIR$ IVDEP
            DO 110 K = ISTAJC( I ), ISTAJC( I + 1 ) - 1
               IG = IGCOLJ( K )
               IF ( IUSED( IG ) == 0 ) THEN
                  WK( IG ) = PI * GRJAC( K )
                  IUSED( IG ) = 1
                  NNONZ2 = NNONZ2 + 1
                  INONZ2( NNONZ2 ) = IG
               ELSE
                  WK( IG ) = WK( IG ) + PI * GRJAC( K )
               END IF
  110       CONTINUE
  120    CONTINUE

!  RESET IUSED TO ZERO.

         DO 130 I = 1, NNONZ2
            IUSED( INONZ2( I ) ) = 0
  130    CONTINUE

!  FORM THE MATRIX-VECTOR PRODUCT Q = A(TRANS) * WK, USING THE
!  ROW-WISE STORAGE OF A.

         NNONNZ = 0
         DO 160 J = 1, NNONZ2
            IG = INONZ2( J )
            IF ( .NOT. GXEQX( IG ) ) THEN

!  IF GROUP IG IS NON TRIVIAL, THERE ARE CONTRIBUTIONS FROM ITS
!  RANK-ONE TERM.

               PI = GSCALE( IG ) * GVALS3( IG ) * WK( IG )
!DIR$ IVDEP
               DO 150 K = ISTAGV( IG ), ISTAGV( IG + 1 ) - 1
                  L = ISVGRP( K )

!  IF Q HAS A NONZERO IN POSITION L, STORE ITS INDEX IN INONNZ.

                  IF ( IUSED( L ) == 0 ) THEN
                     Q( L ) = PI * GRJAC( IVALJR( K ) )
                     IUSED( L ) = 1
                     NNONNZ = NNONNZ + 1
                     INONNZ( NNONNZ ) = L
                  ELSE
                     Q( L ) = Q( L ) + PI * GRJAC( IVALJR( K ) )
                  END IF
  150          CONTINUE
            END IF
  160    CONTINUE
      END IF
      IF ( .NOT. ALLLIN ) THEN

! ======================= SECOND-ORDER TERMS =======================

!  NOW CONSIDER THE PRODUCT OF P WITH THE SECOND ORDER TERMS (THAT IS.
! (2ND DERIVATIVES OF THE ELEMENTS). AGAIN, FOR MAXIMUM EFFICIENCY,
!  THE PRODUCT IS FORMED IN DIFFERENT WAYS IF P IS SPARSE OR DENSE.

! --------------------- CASE 1. P IS NOT SPARSE ---------------------

         IF ( DENSEP ) THEN
            DO 280 IELL = 1, NGEL
               IEL = IELING( IELL )
               IG = ISLGRP( IELL )
               NVAREL = ISTAEV( IEL + 1 ) - ISTAEV( IEL )
               IF ( GXEQX( IG ) ) THEN
                  GI = GSCALE( IG ) * ESCALE( IELL )
               ELSE
                  GI = GSCALE( IG ) * ESCALE( IELL ) * GVALS2( IG ) 
               END IF
               ISWKSP( IEL ) = NBPROD
               IF ( INTREP( IEL ) ) THEN

!  THE IEL-TH ELEMENT HESSIAN HAS AN INTERNAL REPRESENTATION.
!  COPY THE ELEMENTAL VARIABLES INTO WK.

                  NULLWK = .TRUE.
                  LL = ISTAEV( IEL )
!DIR$ IVDEP
                  DO 210 II = 1, NVAREL
                     PI = P( IELVAR( LL ) )
                     WK( II ) = PI
                     IF ( PI /= ZERO ) NULLWK = .FALSE.
                     LL = LL + 1
  210             CONTINUE

!  SKIP THE ELEMENT IF WK IS NULL.

                  IF ( NULLWK ) GO TO 280

!  FIND THE INTERNAL VARIABLES, WKB.

                  NIN = INTVAR( IEL + 1 ) - INTVAR( IEL )
                  CALL RANGE ( IEL, .FALSE., WK, WKB, NVAREL, NIN,  &
                               ITYPEE( IEL ), NVAREL, NIN )

!  MULTIPLY THE INTERNAL VARIABLES BY THE ELEMENT HESSIAN AND PUT THE
!  PRODUCT IN WKC. CONSIDER THE FIRST COLUMN OF THE ELEMENT HESSIAN.

                  IELHST = ISTADH( IEL )
                  PI = GI * WKB( 1 )
!DIR$ IVDEP
                  DO 220 IROW = 1, NIN
                     IJHESS = ISYMMH( 1, IROW ) + IELHST
                     WKC( IROW ) = PI * HUVALS( IJHESS )
  220             CONTINUE

!  NOW CONSIDER THE REMAINING COLUMNS OF THE ELEMENT HESSIAN.

                  DO 240 JCOL = 2, NIN
                     PI = GI * WKB( JCOL )
                     IF ( PI /= ZERO ) THEN
!DIR$ IVDEP
                        DO 230 IROW = 1, NIN
                           IJHESS = ISYMMH( JCOL, IROW ) + IELHST
                           WKC( IROW ) = WKC( IROW ) + PI * &
                                         HUVALS( IJHESS )
  230                   CONTINUE
                     END IF
  240             CONTINUE

!  SCATTER THE PRODUCT BACK ONTO THE ELEMENTAL VARIABLES, WK.

                  CALL RANGE ( IEL, .TRUE., WKC, WK, NVAREL, NIN,  &
                               ITYPEE( IEL ), NIN, NVAREL )

!  ADD THE SCATTERED PRODUCT TO Q.

                  LL = ISTAEV( IEL )
!DIR$ IVDEP
                  DO 250 II = 1, NVAREL
                     L = IELVAR( LL )
                     Q( L ) = Q( L ) + WK( II )
                     LL = LL + 1
  250             CONTINUE
               ELSE

!  THE IEL-TH ELEMENT HESSIAN HAS NO INTERNAL REPRESENTATION.

                  LTHVAR = ISTAEV( IEL ) - 1
                  IELHST = ISTADH( IEL )
                  DO 270 JCOL = 1, NVAREL
                     PI = GI * P( IELVAR( LTHVAR + JCOL ) )
                     IF ( PI /= ZERO ) THEN
!DIR$ IVDEP
                        DO 260 IROW = 1, NVAREL
                           IJHESS = ISYMMH( JCOL, IROW ) + IELHST
                           L = IELVAR( LTHVAR + IROW )
                           Q( L ) = Q( L ) + PI * HUVALS( IJHESS )
  260                   CONTINUE
                     END IF
  270             CONTINUE
               END IF
  280       CONTINUE
         ELSE

! -------------------- CASE 2. P IS SPARSE ------------------------

            DO 380 J = NVAR1, NVAR2

!  CONSIDER EACH NONZERO COMPONENT OF P SEPARATELY.

               I = IVAR( J )
               IPT = ISPTRS( I )
               IF ( IPT >= 0 ) THEN

!  THE INDEX OF THE I-TH COMPONENT LIES IN THE IEL-TH NONLINEAR ELEMENT.

                  IELL = ISELTS( I )
  300             CONTINUE

!  CHECK TO ENSURE THAT THE IEL-TH ELEMENT HAS NOT ALREADY BEEN USED.

                  IF ( ISWKSP( IELL ) < NBPROD ) THEN
                     ISWKSP( IELL ) = NBPROD
                     IEL = IELING( IELL )
                     NVAREL = ISTAEV( IEL + 1 ) - ISTAEV( IEL )
                     IG = ISLGRP( IELL )
                     IF ( GXEQX( IG ) ) THEN
                        GI = GSCALE( IG ) * ESCALE( IELL )
                     ELSE
                        GI = GSCALE( IG ) * ESCALE( IELL ) * GVALS2( IG) 
                     END IF
                     IF ( INTREP( IEL ) ) THEN

!  THE IEL-TH ELEMENT HESSIAN HAS AN INTERNAL REPRESENTATION.
!  COPY THE ELEMENTAL VARIABLES INTO WK.

                        LL = ISTAEV( IEL )
!DIR$ IVDEP
                        DO 310 II = 1, NVAREL
                           WK( II ) = P( IELVAR( LL ) )
                           LL = LL + 1
  310                   CONTINUE

!  FIND THE INTERNAL VARIABLES.

                        NIN = INTVAR( IEL + 1 ) - INTVAR( IEL )
                        CALL RANGE ( IEL, .FALSE., WK, WKB, NVAREL, NIN, &
                                    ITYPEE( IEL ), NVAREL, NIN )

!  MULTIPLY THE INTERNAL VARIABLES BY THE ELEMENT HESSIAN AND PUT THE 
!  PRODUCT IN WKB. CONSIDER THE FIRST COLUMN OF THE ELEMENT HESSIAN.

                        IELHST = ISTADH( IEL )
                        PI = GI * WKB( 1 )
!DIR$ IVDEP
                        DO 320 IROW = 1, NIN
                           IJHESS = ISYMMH( 1, IROW ) + IELHST
                           WKC( IROW ) = PI * HUVALS( IJHESS )
  320                   CONTINUE

!  NOW CONSIDER THE REMAINING COLUMNS OF THE ELEMENT HESSIAN.

                        DO 340 JCOL = 2, NIN
                           PI = GI * WKB( JCOL )
                           IF ( PI /= ZERO ) THEN
!DIR$ IVDEP
                              DO 330 IROW = 1, NIN
                                 IJHESS = ISYMMH( JCOL, IROW ) + IELHST
                                 WKC( IROW ) = WKC( IROW ) + PI * &
                                             HUVALS( IJHESS )
  330                         CONTINUE
                           END IF
  340                   CONTINUE

!  SCATTER THE PRODUCT BACK ONTO THE ELEMENTAL VARIABLES, WK.

                        CALL RANGE ( IEL, .TRUE., WKC, WK, NVAREL, NIN, &
                                     ITYPEE( IEL ), NIN, NVAREL )

!  ADD THE SCATTERED PRODUCT TO Q.

                        LL = ISTAEV( IEL )
!DIR$ IVDEP
                        DO 350 II = 1, NVAREL
                           L = IELVAR( LL )

!  IF Q HAS A NONZERO IN POSITION L, STORE ITS INDEX IN INONNZ.

                           IF ( ABS( WK( II ) ) > TINY ) THEN
                              IF ( IUSED( L ) == 0 ) THEN
                                 Q( L ) = WK( II )
                                 IUSED( L ) = 1
                                 NNONNZ = NNONNZ + 1
                                 INONNZ( NNONNZ ) = L
                              ELSE
                                 Q( L ) = Q( L ) + WK( II )
                              END IF
                           END IF
                           LL = LL + 1
  350                   CONTINUE
                     ELSE

!  THE IEL-TH ELEMENT HESSIAN HAS NO INTERNAL REPRESENTATION.

                        LTHVAR = ISTAEV( IEL ) - 1
                        IELHST = ISTADH( IEL )
                        DO 370 JCOL = 1, NVAREL
                           PI = GI * P( IELVAR( LTHVAR + JCOL ) )
                           IF ( PI /= ZERO ) THEN
!DIR$ IVDEP
                              DO 360 IROW = 1, NVAREL
                                 IJHESS = ISYMMH( JCOL, IROW ) + IELHST

!  IF Q HAS A NONZERO IN POSITION L, STORE ITS INDEX IN INONNZ.

                                 IF ( ABS( HUVALS( IJHESS ) ) &
                                    > TINY ) THEN
                                    L = IELVAR( LTHVAR + IROW )
                                    IF ( IUSED( L ) == 0 ) THEN
                                       Q( L ) = PI * HUVALS( IJHESS )
                                       IUSED( L ) = 1
                                       NNONNZ = NNONNZ + 1
                                       INONNZ( NNONNZ ) = L
                                    ELSE
                                       Q( L ) = Q( L ) + PI * &
                                                HUVALS( IJHESS )
                                    END IF
                                 END IF
  360                         CONTINUE
                           END IF
  370                   CONTINUE
                     END IF
                  END IF

!  CHECK TO SEE IF THERE ARE ANY FURTHER ELEMENTS WHOSE VARIABLES
!  INCLUDE THE I-TH VARIABLE.

                  IF ( IPT > 0 ) THEN
                     IELL = ISELTS( IPT )
                     IPT = ISPTRS( IPT )
                     GO TO 300
                  END IF
               END IF
  380       CONTINUE
         END IF
      END IF

! ==================== THE PRODUCT IS COMPLETE =======================

!  RESET IUSED TO ZERO.

      IF ( .NOT. DENSEP ) THEN
         DO 390 I = 1, NNONNZ
            IUSED( INONNZ( I ) ) = 0
  390    CONTINUE
      END IF
      RETURN

!  END OF HSPRD.

      END
