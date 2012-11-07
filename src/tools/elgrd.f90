! ( Last modified on 23 Dec 2000 at 22:01:38 )

!  ** FOR THE CRAY 2, LINES STARTING CDIR$ IVDEP TELL THE COMPILER TO
!     IGNORE POTENTIAL VECTOR DEPENDENCIES AS THEY ARE KNOWN TO BE O.K.

      SUBROUTINE ELGRD( n, ng, FIRSTG, ICNA, licna, ISTADA, lstada, &
                         IELING, leling, ISTADG, lstadg, ITYPEE, litype, &
                         ISTAEV, lstaev, IELVAR, lelvar, INTVAR, lntvar,  &
                         ISVGRP, lnvgrp, ISTAJC, lnstjc, ISTAGV, lnstgv,  &
                         A, la, GVALS2, lgvals, GUVALS, lguval, GRAD, &
                         GSCALE, lgscal, ESCALE, lescal, GRJAC, lngrjc,  &
                         WKPVAR, WKEVAR, lnwkev, gxeqx, lgxeqx, intrep,  &
                         lintre, RANGE )
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n, ng, la, litype
      INTEGER :: licna, lstada, leling, lstadg, lstaev, lelvar
      INTEGER :: lntvar, lnvgrp, lnstjc, lnstgv, lgvals, lnwkev
      INTEGER :: lgscal, lescal, lgxeqx, lintre, lguval, lngrjc
      LOGICAL :: FIRSTG
      INTEGER :: ICNA( licna ), ISTADA( lstada ), IELING( leling )
      INTEGER :: ISTADG( lstadg ), ISTAEV( lstaev )
      INTEGER :: IELVAR( lelvar ), INTVAR( lntvar )
      INTEGER :: ISVGRP( lnvgrp ), ISTAJC( lnstjc )
      INTEGER :: ISTAGV( lnstgv ), ITYPEE( litype )
      REAL ( KIND = wp ) :: A( la ), GVALS2( lgvals ), GUVALS( lguval ),  &
                       GRAD( n ), GSCALE( lgscal ), ESCALE( lescal ),  &
                       GRJAC( lngrjc ), WKPVAR( n ), WKEVAR( lnwkev )
      LOGICAL :: GXEQX( lgxeqx ), INTREP( lintre )
      EXTERNAL :: RANGE 

!  CALCULATE THE GRADIENT OF EACH GROUP, GRJAC, AND THE GRADIENT OF THE
!  OBJECTIVE FUNCTION, GRAD.

!  NICK GOULD, 4TH JULY 1990.
!  FOR CGT PRODUCTIONS.

!  LOCAL VARIABLES.

      INTEGER :: i, iel, ig, ig1, ii,J, jj, k, l, ll
      INTEGER :: nin, nvarel, nelow, nelup, istrgv, iendgv
      REAL ( KIND = wp ) :: zero, gi, scalee
      LOGICAL :: NONTRV
!D    EXTERNAL         SETVL, SETVI

!  SET CONSTANT REAL PARAMETERS.

      PARAMETER ( zero = 0.0_wp )

!  INITIALIZE THE GRADIENT AS ZERO.

      CALL SETVL( n, GRAD, 1, zero )

!  CONSIDER THE IG-TH GROUP.

      DO 160 ig = 1, ng
         ig1 = ig + 1
         istrgv = ISTAGV( ig )
         iendgv = ISTAGV( ig1 ) - 1
         nelow = ISTADG( ig )
         nelup = ISTADG( ig1 ) - 1
         NONTRV = .NOT. GXEQX( ig )

!  COMPUTE THE FIRST DERIVATIVE OF THE GROUP.

         IF ( NONTRV ) THEN
            gi = GSCALE( ig ) * GVALS2( ig ) 
         ELSE
            gi = GSCALE( ig )
         END IF 

!  THIS is THE FIRST GRADIENT EVALUATION OR THE GROUP HAS NONLINEAR
!  ELEMENTS.

         IF ( FIRSTG .OR. nelow <= nelup ) THEN
      CALL SETVI( iendgv - istrgv + 1, WKPVAR, ISVGRP( istrgv ), &
                         zero )

!  LOOP OVER THE GROUP'S NONLINEAR ELEMENTS.

            DO 30 ii = nelow, nelup
               iel = IELING( ii )
               k = INTVAR( iel )
               l = ISTAEV( iel )
               nvarel = ISTAEV( iel + 1 ) - l
               scalee = ESCALE( ii )
               IF ( INTREP( iel ) ) THEN

!  THE IEL-TH ELEMENT HAS AN INTERNAL REPRESENTATION.

                  nin = INTVAR( iel + 1 ) - k
                  CALL RANGE ( iel, .TRUE., GUVALS( k ), &
                               WKEVAR, nvarel, nin,  &
                               ITYPEE( iel ), nin, nvarel )
!DIR$ IVDEP
                  DO 10 i = 1, nvarel
                     j = IELVAR( l )
                     WKPVAR( j ) = WKPVAR( j ) + scalee * WKEVAR( i )
                     l = l + 1
   10             CONTINUE
               ELSE

!  THE IEL-TH ELEMENT HAS NO INTERNAL REPRESENTATION.

!DIR$ IVDEP
                  DO 20 i = 1, nvarel
                     j = IELVAR( l )
                     WKPVAR( j ) = WKPVAR( j ) + scalee * GUVALS( k )
                     k = k + 1
                     l = l + 1
   20             CONTINUE
               END IF
   30       CONTINUE

!  INCLUDE THE CONTRIBUTION FROM THE LINEAR ELEMENT.

!DIR$ IVDEP
            DO 40 k = ISTADA( ig ), ISTADA( ig1 ) - 1
               j = ICNA( k )
               WKPVAR( j ) = WKPVAR( j ) + A( k )
   40       CONTINUE

!  FIND THE GRADIENT OF THE GROUP.

            IF ( NONTRV ) THEN

!  THE GROUP is NON-TRIVIAL.

!DIR$ IVDEP
               DO 50 i = istrgv, iendgv
                  ll = ISVGRP( i )
                  GRAD( ll ) = GRAD( ll ) + gi * WKPVAR( ll )

!  AS THE GROUP is NON-TRIVIAL, ALSO STORE THE NONZERO ENTRIES OF THE
!  GRADIENT OF THE FUNCTION IN GRJAC.

                  jj = ISTAJC( ll )
                  GRJAC( jj ) = WKPVAR( ll )

!  INCREMENT THE ADDRESS FOR THE next NONZERO IN THE COLUMN OF
!  THE JACOBIAN FOR VARIABLE LL.

                  ISTAJC( ll ) = jj + 1
   50          CONTINUE
            ELSE

!  THE GROUP is TRIVIAL.

!DIR$ IVDEP
               DO 60 i = istrgv, iendgv
                  ll = ISVGRP( i )
                  GRAD( ll ) = GRAD( ll ) + gi * WKPVAR( ll )
   60          CONTINUE
            END IF

!  THIS is NOT THE FIRST GRADIENT EVALUATION AND THERE is ONLY A LINEAR
!  ELEMENT.

         ELSE

!  ADD THE GRADIENT OF THE LINEAR ELEMENT TO THE OVERALL GRADIENT

!DIR$ IVDEP
            DO 130 k = ISTADA( ig ), ISTADA( ig1 ) - 1
               ll = ICNA( k )
               GRAD( ll ) = GRAD( ll ) + gi * A( k )
  130       CONTINUE

!  THE GROUP is NON-TRIVIAL; INCREMENT THE STARTING ADDRESSES FOR
!  THE GROUP  USED BY EACH VARIABLE IN THE (UNCHANGED) LINEAR
!  ELEMENT TO AVOID RESETTING THE NONZEROS IN THE JACOBIAN.

            IF ( NONTRV ) THEN
!DIR$ IVDEP
               DO 140 i = istrgv, iendgv
                  ll = ISVGRP( i )
                  ISTAJC( ll ) = ISTAJC( ll ) + 1
  140          CONTINUE
            END IF
         END IF
  160 CONTINUE

!  RESET THE STARTING ADDRESSES FOR THE LISTS OF GROUP  USING
!  EACH VARIABLE TO THEIR VALUES ON ENTRY.

      DO 170 i = n, 2, - 1
         ISTAJC( i ) = ISTAJC( i - 1 )
  170 CONTINUE
      ISTAJC( 1 ) = 1
      RETURN

!  END OF ELGRD.

      END
