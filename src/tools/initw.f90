! ( Last modified on 23 Dec 2000 at 22:01:38 )

!  ** FOR THE CRAY 2, LINES STARTING 'CDIR$ IVDEP' TELL THE COMPILER TO
!     IGNORE POTENTIAL VECTOR DEPENDENCIES AS THEY ARE KNOWN TO BE O.K.

      SUBROUTINE INITW( n, ng, nel, ieling, leling, istadg, lstadg, &
                         ielvar, lelvar, istaev, lstaev, intvar, lntvar, &
                         istadh, lstadh, icna, licna, istada, lstada, &
                         itypee, litype, &
                         gxeqx, lgxeqx, intrep, lintre, lfuval, ALTRIV, &
                         DIRECT, FDGRAD, lfxi, lgxi, lhxi, lggfx, &
                         ldx, lgrjac, lqgrad, lbreak, lp, lxcp, &
                         lx0, lgx0, ldeltx, lbnd, lwkstr, lsptrs, &
                         lselts, lindex, lswksp, lstagv, lstajc, liused, &
                         lfreec, lnnonz, lnonz2, lsymmd, lsymmh, lslgrp, &
                         lsvgrp, lgcolj, lvaljr, lsend, lnptrs, lnelts, &
                         lnndex, lnwksp, lnstgv, lnstjc, lniuse, lnfrec, &
                         lnnnon, lnnno2, lnsymd, lnsymh, lnlgrp, lnvgrp, &
                         lngclj, lnvljr, lnqgrd, lnbrak, lnp, lnbnd, &
                         lnfxi, lngxi, lnguvl, lnhxi, lnhuvl, lnggfx, &
                         lndx, lngrjc, liwk2, lwk2, maxsin, ninvar, &
                         ntype, nsets, maxsel, lstype, lsswtr, lssiwt, &
                         lsiwtr, lswtra, lntype, lnswtr, lnsiwt, lniwtr, &
                         lnwtra, lsiset, lssvse, lniset, lnsvse, RANGE, &
                         IWK, liwk, WK,  lwk, iprint, iout, INFORM)
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

! --------------------------------------------------------------------

!  COMPUTE THE STARTING ADDRESSES FOR THE PARTITIONS OF THE WORKSPACE
!  ARRAYS FUVALS, IWK AND WK. ALSO FILL RELEVANT PORTIONS OF IWK.
!  THE ADDRESSES FOR IWK AND WK ARE DESCRIBED HERE; THOSE FOR FUVALS
!  ARE DESCRIBED IN THE INTRODUCTORY COMMENTS TO SUBROUTINE SBMIN.

!  NICK GOULD, 20TH JUNE 1990.
!  FOR CGT PRODUCTIONS.

! --------------------------------------------------------------------

      INTEGER :: n, ng, nel, inform, litype, lwk, liwk
      INTEGER :: lfuval, lelvar, lstaev, lstadh, leling, lntvar
      INTEGER :: lgxeqx, lstadg, licna, lstada, iout, iprint
      INTEGER :: lfxi, lgxi, lhxi, lggfx, ldx, lgrjac
      INTEGER :: lqgrad, lbreak, lp, lxcp, lbnd, lwkstr
      INTEGER :: lgx0, ldeltx, lsptrs, lselts, lindex, lx0
      INTEGER :: lswksp, lstagv, lstajc, liused, lfreec, lintre
      INTEGER :: lnnonz, lnonz2, lsymmd, lsymmh, lsiset, lssvse
      INTEGER :: lslgrp, lsvgrp, lgcolj, lvaljr, lsend
      INTEGER :: lnptrs, lnelts, lnndex, lnwksp, lnstgv, lnstjc
      INTEGER :: lniuse, lnnnon, lnnno2, lnsymd, lnsymh, lnlgrp
      INTEGER :: lnvgrp, lngclj, lnvljr, lnfrec, lniset, lnsvse
      INTEGER :: lnqgrd, lnbrak, lnp, lnbnd,  lnfxi, lngxi
      INTEGER :: lnguvl, lnhxi, lnhuvl, lnggfx, lndx, lngrjc
      INTEGER :: liwk2, lwk2, ninvar, maxsel, maxsin, ntype
      INTEGER :: lstype, lsswtr, lssiwt, lsiwtr, lswtra, nsets
      INTEGER :: lntype, lnswtr, lnsiwt, lniwtr, lnwtra
      EXTERNAL :: RANGE 
      LOGICAL :: DIRECT, ALTRIV, FDGRAD
      INTEGER :: IELVAR( lelvar ), ISTAEV( lstaev )
      INTEGER :: ISTADH( lstadh ), IWK ( liwk )
      INTEGER :: INTVAR( lntvar ), ISTADG( lstadg )
      INTEGER :: icna ( licna ), ISTADA( lstada )
      INTEGER :: IELING( leling ), ITYPEE( litype )
      REAL ( KIND = wp ) :: WK ( lwk )
      LOGICAL :: gxeqx ( lgxeqx ), INTREP( lintre )

!  LOCAL VARIABLES.

      INTEGER :: i, j, k, l, nvargp, iielts, ientry, ig, is
      INTEGER :: lend, nsizeh, nel1, ngel, ng1, iel
      INTEGER :: lwfree, liwfre, nelvr, lw1, liwfro, lwfreo
      INTEGER :: iell, itype, isofar, istarj, ivarp1, ivar
      INTEGER :: jset, inext, newvar, newset, ipt, istrt
      INTEGER :: ninvr, ii, jj, kk, ll
      LOGICAL :: NONTRV, ALLLIN, VRUSED
      REAL ( KIND = wp ) :: zero, one
      PARAMETER ( zero = 0.0_wp, one = 1.0_wp )

!  EXTERNAL SUBROUTINES AND FUNCTIONS USED.

!D    EXTERNAL         SYMMH

!  INTRINSIC FUNCTIONS.

      INTRINSIC        MAX

!  SET CONSTANTS.

      nel1 = nel + 1
      ng1 = ng + 1
      ngel = ISTADG( ng1 ) - 1
      ALLLIN = nel == 0

!  SET UP intvar, THE STARTING ADDRESSES FOR THE ELEMENT GRADIENTS
!  WITH RESPECT TO THEIR INTERNAL VARIABLES. ALSO COMPUTE maxsin,
!  THE MAXIMUM NUMBER OF INTERNAL VARIABLES IN AN ELEMENT.

      IF ( .NOT. ALLLIN ) THEN
         k = INTVAR( 1 )
         maxsin = k
         INTVAR( 1 ) = nel1
         DO 10 iel = 2, nel
            l = INTVAR( iel )
            INTVAR( iel ) = INTVAR( iel - 1 ) + k
            k = l
            maxsin = MAX( maxsin, k )
   10    CONTINUE
         INTVAR( nel1 ) = INTVAR( nel ) + k
      ELSE
         INTVAR( 1 ) = 1
         maxsin = 0
      END IF

!  COMPUTE THE TOTAL NUMBER OF INTERNAL VARIABLES.

      ninvar = INTVAR( nel1 ) - INTVAR( 1 )

!  CALCULATE THE LENGTH, iielts, OF WORKSPACE REQUIRED TO
!  DETERMINE WHICH ELEMENTS USE EACH OF THE VARIABLES.
!  ALSO FIND THE MAXIMUM NUMBER OF VARIABLES IN AN ELEMENT, MAXSEL.
!  THIS is A DUMMY RUN FOR LOOP 130 MERELY TO CALCULATE THE
!  SPACE REQUIRED.

      IF ( liwk < n ) THEN
         WRITE( iout, 2030 ) n - liwk
         inform = 4
         RETURN
      END IF

!  IWK( i ) WILL BE USED AS A LIST OF LINKS CHAINING THE ELEMENTS USING
!  VARIABLE I. IF IWK( i ) is NEGATIVE, THE LIST is EMPTY.

!DIR$ IVDEP
      DO 20 i = 1, n
         IWK( i ) = - 1
   20 CONTINUE
      iielts = n
      maxsel = 0
      IF ( .NOT. ALLLIN ) THEN

!  LOOP OVER THE GROUP, CONSIDERING EACH NONLINEAR ELEMENT IN TURN.

         DO 50 i = 1, ngel
            iel = IELING( i )
            maxsel = MAX( maxsel, ISTAEV( iel + 1 ) - ISTAEV( iel ) )

!  LOOP ON THE VARIABLES FROM THE I-TH ELEMENT.

            DO 40 k = ISTAEV( iel ), ISTAEV( iel + 1 ) - 1
               ientry = IELVAR( k )
               IF ( IWK( ientry ) >= 0 ) THEN

!  IF WE HAVE REACHED THE END OF THE LIST OF THE ELEMENTS USING
!  THE VARIABLE IELVAR( k ), ADD THE IEL-TH ELEMENT TO IT.
!  OTHERWISE, FIND THE next ENTRY IN THE LIST.

   30             CONTINUE
                  IF ( IWK( ientry ) > 0 ) THEN
                     ientry = IWK( ientry )
                     GO TO 30
                  ELSE
                     iielts = iielts + 1
                     IF ( iielts > liwk ) THEN
                        WRITE( iout, 2030 ) iielts - liwk
                        inform = 4
                        RETURN
                     END IF
                     IWK( ientry ) = iielts
                     IWK( iielts ) = 0
                  END IF
               ELSE

!  THE LIST OF ELEMENTS INVOLVING THE VARIABLE IELVAR( k ) WAS
!  PREVIOUSLY EMPTY. INDICATE THAT THE LIST HAS NOW BEEN STARTED AND
!  THAT ITS END HAS BEEN REACHED.

                  IWK( ientry ) = 0
               END IF
   40       CONTINUE
   50    CONTINUE
      END IF

! -- CALCULATE THE STARTING ADDRESSES FOR THE INTEGER WORKSPACE. --

!  IWK( lsptrs + j ), j = 1, ..., iielts,  WILL CONTAIN THE LINKS FOR
!  THE LISTS OF NONLINEAR ELEMENTS WHICH USE EACH VARIABLE.

      lsptrs = 0

!  IWK( lselts + j ), j = 1, ..., iielts, WILL CONTAIN THE LISTS OF
!  NONLINEAR ELEMENTS CORRESPONDING TO THE PREVIOUSLY MENTIONED LINKS.

      lselts = lsptrs + iielts

!  IWK( lindex + j ), j = 1, ..., n, WILL CONTAIN THE STATUS OF THE
!  J-TH VARIABLE AS THE CURRENT ITERATION PROGRESSES. POSSIBLE VALUES
!  ARE 0 IF THE VARIABLE LIES AWAY FROM ITS BOUNDS, 1 AND 2 IF IT LIES
!  ON ITS LOWER OR UPPER BOUNDS (RESPECTIVELY) - THESE MAY BE PROBLEM
!  BOUNDS OR TRUST REGION BOUNDS, AND 3 IF THE VARIABLE is FIXED.

      lindex = lselts + iielts

!  IWK( lswksp + j ), j = 1, ..., MAX( ngel, nel, n + n ), is USED FOR
!  WORKSPACE BY THE MATRIX-VECTOR PRODUCT SUBROUTINE HESPRD.

      lswksp = lindex + n

!  IWK( liused + j ), j = 1, ..., MAX( n, ng ) WILL BE USED AS
!  WORKSPACE BY THE MATRIX-VECTOR PRODUCT SUBROUTINE HESPRD.

      IF ( DIRECT ) THEN
         liused = lswksp + MAX( ngel, nel, n + n )
      ELSE
         liused = lswksp + MAX( ngel, nel, n )
      END IF

!  IWK( LFREED + j ), j = 1, ..., NFREEC WILL GIVE THE INDICES OF THE
!  VARIABLES WHICH ARE CONSIDERED TO BE FREE FROM THEIR BOUNDS AT THE
!  CURRENT GENERALIZED CAUCHY POINT.

      lfreec = liused + MAX( n, ng )

!  IWK( lnnonz + j ), j = 1, ..., NNNONZ WILL GIVE THE INDICES OF THE
!  NONZEROS IN THE VECTOR OBTAINED AS A RESULT OF THE MATRIX-VECTOR
!  PRODUCT FROM SUBROUTINE HESPRD.

      lnnonz = lfreec + n

!  IWK( lnonz2 + j ), j = 1, ..., ng, WILL BE USED AS FURTHER
!  WORKSPACE BY THE MATRIX-VECTOR PRODUCT SUBROUTINE HESPRD.

      lnonz2 = lnnonz + n

!  IWK( lsymmd + j ), j = 1, ..., maxsin, WILL GIVE THE LOCATION OF
!  THE J-TH DIAGONAL OF A maxsin BY maxsin SYMMETRIC MATRIX IN AN
!  UPPER TRIANGULAR STORAGE SCHEME.

      lsymmd = lnonz2 + ng

!  IWK( lsymmh + ( i - 1 ) * maxsin + j ), i, j = 1, ..., maxsin,
!  WILL GIVE THE LOCATION OF THE (I,J)-TH ENTRY OF A maxsin BY
!  maxsin SYMMETRIC MATRIX IN AN UPPER TRIANGULAR STORAGE SCHEME.

      lsymmh = lsymmd + maxsin

!  IWK( lslgrp + j ), j = 1, ..., ngel, WILL CONTAIN THE NUMBER OF
!  THE GROUP WHICH USES NONLINEAR ELEMENT J.

      lslgrp = lsymmh + maxsin * maxsin

!  IWK( lstajc + j ), j = 1, ..., n, WILL CONTAIN THE STARTING
!  ADDRESSES FOR THE LIST OF NONTRIVIAL GROUP  WHICH USE THE
!  J-TH VARIABLE. IWK( lstajc + n + 1 ) WILL POINT TO THE FIRST FREE
!  LOCATION IN IWK AFTER THE LIST OF NONTRIVIAL GROUP  FOR THE
!  N-TH VARIABLE.

      lstajc = lslgrp + ngel

!  IWK( lstagv + j ), j = 1, ..., ng, WILL CONTAIN THE STARTING
!  ADDRESSES FOR THE LIST OF VARIABLES WHICH OCCUR IN THE J-TH GROUP.
!  IWK( lstagv + ng + 1 ) WILL POINT TO THE FIRST FREE LOCATION
!  IN IWK AFTER THE LIST OF VARIABLES FOR THE NG-TH GROUP.

      lstagv = lstajc + n + 1

!  IWK( lsvgrp + j ), j = 1, ..., nvargp, WILL CONTAIN THE INDICES OF
!  THE VARIABLES WHICH ARE USED BY EACH GROUP IN TURN. THOSE FOR GROUP i
!  OCCUR IN LOCATIONS IWK( lstagv + i ) TO IWK( lstagv + i + 1 ) - 1.

      lsvgrp = lstagv + ng + 1

!  CHECK THAT THERE is SUFFICIENT WORKSPACE.

      IF ( lsvgrp > liwk ) THEN
         WRITE( iout, 2030 ) lsvgrp - liwk
         inform = 4
         RETURN
      END IF

!  IWK( lgcolj + j ), j = 1, ..., nvargp, WILL CONTAIN THE INDICES OF
!  THE NONTRIVIAL GROUP  WHICH USE EACH VARIABLE IN TURN. THOSE FOR
!  VARIABLE i OCCUR IN LOCATIONS IWK( lstajc + i ) TO
!  IWK( lstajc + i + 1 ) - 1.

      lgcolj = lsvgrp + 1

!  DETERMINE WHICH ELEMENTS USE EACH VARIABLE. INITIALIZATION.

      IF ( .NOT. ALLLIN ) THEN

!  IWK( i ) WILL BE USED AS A LIST OF LINKS CHAINING THE ELEMENTS USING
!  VARIABLE I. IF IWK( i ) is NEGATIVE, THE LIST is EMPTY.

!DIR$ IVDEP
         DO 100 i = 1, n
            IWK( i ) = - 1
  100    CONTINUE
         iielts = n

!  LOOP OVER THE GROUP, CONSIDERING EACH NONLINEAR ELEMENT IN TURN.

         DO 130 i = 1, ngel
            iel = IELING( i )

!  LOOP ON THE VARIABLES OF THE I-TH ELEMENT.

            DO 120 k = ISTAEV( iel ), ISTAEV( iel + 1 ) - 1
               ientry = IELVAR( k )
               IF ( IWK( ientry ) >= 0 ) THEN

!  IF WE HAVE REACHED THE END OF THE LIST OF THE ELEMENTS USING
!  THE VARIABLE IELVAR( k ), ADD THE I-TH ELEMENT TO IT AND
!  RECORD THAT THE END OF THE LIST HAS OCCURED.
!  OTHERWISE, FIND THE next ENTRY IN THE LIST.

  110             CONTINUE
                  IF ( IWK( ientry ) > 0 ) THEN
                     ientry = IWK( ientry )
                     GO TO 110
                  ELSE
                     iielts = iielts + 1
                     IWK( ientry ) = iielts
                     IWK( lselts + iielts ) = i
!                    IWK( lselts + iielts ) = iel
                     IWK( iielts ) = 0
                  END IF
               ELSE

!  THE LIST OF ELEMENTS INVOLVING THE VARIABLE IELVAR( k ) WAS
!  PREVIOUSLY EMPTY. INDICATE THAT THE LIST HAS NOW BEEN STARTED,
!  RECORD THE ELEMENT WHICH CONTAINS THE VARIABLE AND INDICATE
!  THAT THE END OF THE LIST HAS BEEN REACHED.

                  IWK( lselts + ientry ) = i
!                 IWK( lselts + ientry ) = iel
                  IWK( ientry ) = 0
               END IF
  120       CONTINUE
  130    CONTINUE
      END IF

!  SET UP SYMMETRIC ADDRESSES FOR THE UPPER TRIANGULAR STORAGE
!  SCHEMES FOR THE ELEMENT HESSIANS.

      IF ( maxsin > 0 ) CALL SYMMH ( maxsin, IWK( lsymmh + 1 ), &
                                        IWK( lsymmd + 1 ) )

!  SET UP THE STARTING ADDRESSES FOR THE ELEMENT HESSIANS
!  WITH RESPECT TO THEIR INTERNAL VARIABLES AND A POINTER BEYOND
!  THE END OF THE SPACE REQUIRED FOR THE HESSIANS.

      lggfx = INTVAR( nel1 )
      IF ( .NOT. ALLLIN ) THEN
         DO 140 i = 1, nel
            ISTADH( i ) = lggfx
            nsizeh = INTVAR( i + 1 ) - INTVAR( i )
            lggfx = lggfx + nsizeh * ( nsizeh + 1 ) / 2
  140    CONTINUE
      END IF
      ISTADH( nel1 ) = lggfx

! -- CALCULATE THE STARTING ADDRESSES FOR THE REAL WORKSPACE. --

!  WK( lqgrad + j ), j = 1, ..., n, WILL CONTAIN THE GRADIENT OF
!  THE QUADRATIC MODEL AT THE CURRENT ESTIMATE OF THE MINIMIZER.

      lqgrad = MAX( MAX( ng, maxsel ) + 2 * maxsin, &
                    n + ninvar + MAX( maxsel, ninvar ) )

!  WK( lbreak + j ), j = 1, ..., n, WILL CONTAIN THE BREAKPOINTS
!  ALONG THE CAUCHY ARC FROM THE CURRENT ESTIMATE OF THE MINIMIZER.

      lbreak = lqgrad + n

!  WK( lp + j ), j = 1, ..., n, WILL CONTAIN THE VECTOR REQUIRED
!  BY THE MATRIX-VECTOR PRODUCT SUBROUTINE HESPRD.

      lp = lbreak + n

!  WK( lxcp + j ), j = 1, ..., n, WILL CONTAIN THE CURRENT
! (APPROXIMATE) CAUCHY POINT.

      lxcp = lp + n

!  WK( lx0 + j ), j = 1, ..., n, WILL CONTAIN THE CURRENT
!  START OF THE CAUCHY SEARCH. THIS FEATURE is ONLY USED IF
!  MORE THAN one CYCLE is USED TO SOLVE THE BQP ACCURATELY.

      lx0 = lxcp + n

!  WK( lgx0 + j ), j = 1, ..., n, WILL CONTAIN THE GRADIENT AT THE
!  CURRENT START OF THE CAUCHY SEARCH. THIS FEATURE is ONLY USED IF
!  MORE THAN one CYCLE is USED TO SOLVE THE BQP ACCURATELY.

      lgx0 = lx0 + n

!  WK( ldeltx + j ), j = 1, ..., n, WILL CONTAIN THE STEP TAKEN
!  DURING A SINGLE CYCLE, IF MORE THAN one CYCLE is USED TO SOLVE
!  THE BQP ACCURATELY.

      ldeltx = lgx0 + n

!  WK( lbnd + 2 * ( j - 1 ) + i ), j = 1, ..., n, i = 1, 2, WILL
!  CONTAIN THE CURRENT LOWER (I=1) AND UPPER (I=2) BOUNDS ON THE
!  VARIABLES DEFINED BY THE INTERSECTION OF THE TRUST REGION WITH
!  THE FEASIBLE BOX.

      lbnd = ldeltx + n

!  WK( lwkstr + j ), j = 1, ..., lwk2, is THE REMAINING REAL
!  WORKSPACE WHICH is FREE FOR OTHER PURPOSES, SUCH AS FOR FORMING
!  THE FACTORIZATION OF THE MODEL HESSIAN, IF REQUIRED.

      lwkstr = lbnd + n + n
      lwk2 = lwk - lwkstr

!  CHECK THAT THERE is SUFFICIENT REAL WORKSPACE.

      IF ( lwkstr > lwk ) THEN
         inform = 5
         WRITE( iout, 2040 ) lwkstr - lwk
         RETURN
      END IF

!  SET THE LENGTH OF EACH PARTITION OF THE REAL WORKSPACE ARRAY
!  FUVALS FOR ARRAY BOUND CHECKING IN CALLS TO OTHER SUBPROGRAMS.

      lnqgrd = lbreak - lqgrad
      lnbrak = lp - lbreak
      lnp = lxcp - lp
      lnbnd = lwkstr - lbnd

!  STORE THE INDICES OF VARIABLES WHICH APPEARS IN EACH GROUP
!  AND HOW MANY GROUP  USE EACH VARIABLE. START BY INITIALIZING
!  COUNTING ARRAYS TO ZERO.

!DIR$ IVDEP
      DO 150 j = 1, n
         IWK( lswksp + j ) = 0
         IWK( lstajc + j + 1 ) = 0
  150 CONTINUE

!  ALTRIV SPECIFIES WHETHER ALL THE GROUP  ARE TRIVIAL.

      ALTRIV = .TRUE.

!  COUNT THE TOTAL NUMBER OF VARIABLES IN ALL THE GROUP, NVARGP.

      nvargp = 0
      IWK( lstagv + 1 ) = 1

!  LOOP OVER THE GROUP . SEE IF THE IG-TH GROUP is TRIVIAL.

      DO 200 ig = 1, ng
         NONTRV = .NOT. GXEQX( ig )

!  CHECK TO SEE IF ALL OF THE GROUP  ARE TRIVIAL.

         IF ( NONTRV ) ALTRIV = .FALSE.

!  LOOP OVER THE NONLINEAR ELEMENTS FROM THE IG-TH GROUP.

         DO 170 k = ISTADG( ig ), ISTADG( ig + 1 ) - 1
            iel = IELING( k )

!  RUN THROUGH ALL THE ELEMENTAL VARIABLES CHANGING THE I-TH ENTRY OF
!  IWK( lswksp ) FROM zero TO one IF VARIABLE i APPEARS IN AN ELEMENT.

            DO 160 j = ISTAEV( iel ), ISTAEV( iel + 1 ) - 1
               i = IELVAR( j )
               IF ( IWK( lswksp + i ) == 0 ) THEN
                  IWK( lswksp + i ) = 1

!  IF THERE is SUFFICIENT ROOM, RECORD THE NONLINEAR VARIABLES FROM
!  THE IG-TH GROUP.

                  IF ( lgcolj > liwk ) THEN
                     WRITE( iout, 2030 ) lgcolj - liwk
                     inform = 4
                     RETURN
                  END IF
                  IWK( lgcolj ) = i
                  lgcolj = lgcolj + 1
                  nvargp = nvargp + 1
               END IF
  160       CONTINUE

!  RECORD THAT NONLINEAR ELEMENT k OCCURS IN GROUP IELGRP(IEL).

            IWK( lslgrp + k ) = ig
  170    CONTINUE

!  CONSIDER VARIABLES WHICH ARISE FROM THE LINEAR ELEMENT.

         DO 180 j = ISTADA( ig ), ISTADA( ig + 1 ) - 1
            i = ICNA( j )
            IF ( IWK( lswksp + i ) == 0 ) THEN
               IWK( lswksp + i ) = 1

!  IF THERE is SUFFICIENT ROOM, RECORD THE LINEAR VARIABLES FROM
!  THE IG-TH GROUP.

               IF ( lgcolj > liwk ) THEN
                  WRITE( iout, 2030 ) lgcolj - liwk
                  inform = 4
                  RETURN
               END IF
               IWK( lgcolj ) = i
               lgcolj = lgcolj + 1
               nvargp = nvargp + 1
            END IF
  180    CONTINUE

!  RESET THE STATUS ARRAY IWK( lswksp ) TO ZERO.

         DO 190 j = IWK( lstagv + ig ), nvargp
            l = IWK( lsvgrp + j )
            IWK( lswksp + l ) = 0

!  RECORD THAT one FURTHER NONTRIVIAL GROUP USES VARIABLE L.

            IF ( NONTRV ) IWK( lstajc + l + 1 ) = &
                          IWK( lstajc + l + 1 ) + 1
  190    CONTINUE

!  RECORD THE STARTING ADDRESS OF THE VARIABLES IN THE next GROUP.

         IWK( lstagv + ig + 1 ) = nvargp + 1
  200 CONTINUE

! -- CONTINUE SETTING STARTING ADDRESSES FOR PARTITIONS OF IWK. --

!  IWK( lvaljr + j ), j = 1, ..., nvargp, WILL CONTAIN THE POSITIONS IN
!  FUVALS (RELATIVE TO THE STARTING ADDRESS LGRJAC) OF THE NONZEROS OF
!  THE JACOBIAN OF THE GROUP  CORRESPONDING TO THE VARIABLES AS ORDERED
!  IN IWK( lsvgrp + j ).

      lvaljr = lgcolj + nvargp

!  lsend GIVES THE TOTAL FIXED AMOUNT OF INTEGER WORKSPACE USED. liwk2
!  GIVES THE AMOUNT OF WORKSPACE WHICH is FREE FOR OTHER PURPOSES, SUCH
!  AS FOR FORMING THE FACTORIZATION OF THE MODEL HESSIAN, IF REQUIRED.

      lsend = lvaljr + nvargp

!  CHECK THAT THERE is SUFFICIENT INTEGER WORKSPACE.

      IF ( lsend > liwk ) THEN
         WRITE( iout, 2030 ) lsend - liwk
         inform = 4
         RETURN
      END IF

!  SET THE LENGTH OF EACH PARTITION OF THE INTEGER WORKSPACE FOR
!  ARRAY BOUND CHECKING IN CALLS TO OTHER SUBPROGRAMS.

      lnptrs = MAX( 1, lselts - lsptrs )
      lnelts = MAX( 1, lindex - lselts )
      lnndex = MAX( 1, lswksp - lindex )
      lnwksp = MAX( 1, liused - lswksp )
      lniuse = MAX( 1, lfreec - liused )
      lnfrec = MAX( 1, lnnonz - lfreec )
      lnnnon = MAX( 1, lnonz2 - lnnonz )
      lnnno2 = MAX( 1, lsymmd - lnonz2 )
      lnsymd = MAX( 1, lsymmh - lsymmd )
      lnsymh = MAX( 1, lslgrp - lsymmh )
      lnlgrp = MAX( 1, lstajc - lslgrp )
      lnstjc = MAX( 1, lstagv - lstajc )
      lnstgv = MAX( 1, lsvgrp - lstagv )
      lnvgrp = MAX( 1, lgcolj - lsvgrp )
      lngclj = MAX( 1, lvaljr - lgcolj )
      lnvljr = MAX( 1, lsend - lvaljr )

!  SET THE STARTING ADDRESSES FOR THE LISTS OF NONTRIVIAL GROUP 
!  WHICH USE EACH VARIABLE IN TURN.

      k = lstajc + 1
      IWK( k ) = 1
      DO 210 i = 2, n + 1
         k = k + 1
         IWK( k ) = IWK( k ) + IWK( k - 1 )
  210 CONTINUE

!  CONSIDER THE IG-TH GROUP IN ORDER TO ASSOCIATE VARIABLES WITH GROUP .

      DO 230 ig = 1, ng
         IF ( .NOT. GXEQX( ig ) ) THEN
            DO 220 i = IWK( lstagv + ig ), IWK( lstagv + ig + 1 ) - 1
                l = lstajc + IWK( lsvgrp + i )

!  RECORD THAT GROUP ig USES VARIABLE IWK( lsvgrp + i ).

               j = IWK( l )
               IWK( lgcolj + j ) = ig

!  STORE THE LOCATIONS IN THE JACOBIAN OF THE GROUP  OF THE NONZEROS
!  CORRESPONDING TO EACH VARIABLE IN THE IG-TH GROUP. INCREMENT THE
!  STARTING ADDRESS FOR THE POINTER TO THE next GROUP USING VARIABLE
!  IWK( lsvgrp + i ).

               IWK( lvaljr + i ) = j
               IWK( l ) = j + 1
  220      CONTINUE
         END IF
  230 CONTINUE

!  RESET THE STARTING ADDRESSES FOR THE LISTS OF GROUP  USING
!  EACH VARIABLE.

      DO 240 i = n, 2, - 1
         l = lstajc + i
         IWK( l ) = IWK( l - 1 )
  240 CONTINUE
      IWK( lstajc + 1 ) = 1

!  INITIALIZE WORKSPACE VALUES FOR SUBROUTINE HESPRD.

!DIR$ IVDEP
      DO 250 j = 1, MAX( n, ng )
         IWK( liused + j ) = 0
  250 CONTINUE

!  DEFINE FURTHER PARTITIONS OF THE WORKSPACE WHENEVER FINITE-
! -DIFFERENCE GRADIENTS ARE USED.

      IF (  FDGRAD ) THEN

! -- CONTINUE SETTING STARTING ADDRESSES FOR PARTITIONS OF IWK. --

!  THE RANGE TRANSFORMATION FOR EACH NONLINEAR ELEMENT is OF A GIVEN
!  TYPE. SUPPOSE THERE ARE ntype NON-TRIVIAL TYPES. IWK( lstype + i )
!  GIVES THE TYPE OF NONLINEAR ELEMENT i FOR i = 1, ...., NEL.

         lstype = lsend

!  THE RANGE TRANSFORMATION FROM ELEMENTAL TO INTERNAL VARIABLES is
!  DEFINED BY A MATRIX W. FOR EACH NON-TRIVIAL TRANSFORMATION, THE
!  MATRIX W is RECORDED. THE INFORMATION FOR THE I-TH TYPE STARTS IN
!  LOCATION lswtra + IWK( lsswtr + i ), i = 1, ...., NTYPE.

         lsswtr = lstype + nel

!  FOR EACH TYPE OF NONLINEAR ELEMENT USING A NONTRIVIAL RANGE
!  TRANSFORMATION, INTEGER INFORMATION ISD ALSO RECORDED.
!  THE INFORMATION FOR THE I-TH TYPE STARTS IN LOCATION
!  lsiwtr + IWK( lssiwt + i ), i = 1, ...., NTYPE.

         lssiwt = lsswtr + nel

!  THE FOLLOWING PIECES OF INTEGER INFORMATION ARE RECORDED ABOUT
!  THE I-TH TYPE OF NONLINEAR ELEMENT:

!    IWK( lsiwtr + IWK( lssiwt + i ) + 1 ):
!            THE NUMBER OF INTERNAL VARIABLES, NINVR.
!    IWK( lsiwtr + IWK( lssiwt + i ) + 2 ):
!            THE NUMBER OF ELEMENTAL VARIABLES, NELVR.
!    IWK( lsiwtr + IWK( lssiwt + i ) + 2 + i ),
!         i = 1, ..., nelvr + NINVR:
!            PIVOT SEQUENCES FOR THE lu FACTORS OF W.

!  AFTER THE FACTORIZATION AND COMPRESSION, ONLY ninvr LINEARLY
!  INDEPENDENT COLUMNS OF W ARE STORED.

         lsiwtr = lssiwt + nel
         IF ( liwk < lsiwtr ) THEN
            WRITE( iout, 2030 ) lsiwtr - liwk
            inform = 4
            RETURN
         END IF

! -- CONTINUE SETTING STARTING ADDRESSES FOR PARTITIONS OF WK. --

!  THE FOLLOWING PIECES OF INTEGER INFORMATION ARE RECORDED ABOUT
!  THE I-TH TYPE OF NONLINEAR ELEMENT:

!    IWK( lswtra + IWK( lsswtr + i ) + i ),
!         i = 1, ..., nelvr * NINVR:  THE MATRIX W STORED BY COLUMNS.

!  AFTER THE FACTORIZATION AND COMPRESSION, ONLY ninvr LINEARLY
!  INDEPENDENT COLUMNS OF W ARE STORED.

         lswtra = lwkstr

! ---------------------------------------------------------------------
!  CONSIDER ONLY ELEMENTS WHICH USE INTERNAL VARIABLES.
! ---------------------------------------------------------------------

         ntype = 0
         lwfree = 1
         liwfre = 1

!  LOOP OVER ALL NONLINEAR ELEMENTS.

         DO 350 iel = 1, nel
            IF ( INTREP( iel ) ) THEN

!  CALCULATE THE RANGE TRANSFORMATION MATRIX W.

               is = ISTAEV( iel )
               ninvr = INTVAR( iel + 1 ) - INTVAR( iel )
               nelvr = ISTAEV( iel + 1 ) - is
               lw1 = lwfree + ninvr * nelvr - 1
               l = lswtra + lw1
!DIR$ IVDEP
               DO 280 i = 1, nelvr
                  WK( l + i ) = zero
  280          CONTINUE
               k = lswtra + lwfree
               is = is - 1
               DO 320 i = 1, nelvr
                  WK( l + i ) = one
                  CALL RANGE ( iel, .FALSE., WK( l + 1 ), &
                               WK( k ), nelvr, ninvr,  &
                               ITYPEE( iel ), nelvr, ninvr )
                  WK( l + i ) = zero
                  k = k + ninvr

!  CHECK TO SEE IF ANY OF THE COLUMNS BELONG TO DUPLICATED VARIABLES.

                  ii = IELVAR( is + i )
                  DO 290 j = 1, i - 1
                     IF ( IELVAR( is + j ) == ii ) GO TO 300
  290             CONTINUE
                  GO TO 320

!  AMALGAMATE COLUMNS FROM DUPLICATE VARIABLES.

  300             CONTINUE
                  kk = lswtra + lwfree + ( j - 1 ) * ninvr - 1
                  ll = k - ninvr - 1
                  DO 310 jj = 1, ninvr
                     WK( kk + jj ) = WK( kk + jj ) + WK( ll + jj )
                     WK( ll + jj ) = zero
  310             CONTINUE
  320          CONTINUE

!  COMPARE THIS TRANSFORMATION MATRIX WITH PREVIOUS ONES.

               DO 340 i = 1, ntype
                  IF ( IWK( lsiwtr + IWK( lssiwt + i ) ) /= ninvr &
                       .OR. IWK( lsiwtr + IWK( lssiwt + i ) + 1 ) &
                       /= nelvr ) GO TO 340
                  DO 330 j = 0, ninvr * nelvr - 1
                     IF ( WK( lswtra + lwfree + j ) /= WK( &
                          lswtra + IWK( lsswtr + i ) + j ) ) GO TO 340
  330             CONTINUE

!  THE TRANSFORMATION is AN EXISTING ONE. RECORD WHICH ONE.

                  IWK( lstype + iel ) = i
                  GO TO 350
  340          CONTINUE

!  ENSURE THAT THERE is SUFFICIENT ROOM.

               IF ( liwfre + 2 + ninvr + nelvr > liwk ) THEN
                  inform = 4
                  WRITE( iout, 2030 ) liwfre + 2 + ninvr + nelvr - liwk
                  RETURN
               END IF
               IF ( lwfree + ninvr * nelvr > lwk ) THEN
                  inform = 5
                  WRITE( iout, 2040 ) lwfree + ninvr * nelvr - lwk
                  RETURN
               END IF

!  THE TRANSFORMATION DEFINES A NEW TYPE. RECORD ITS DETAILS.

               ntype = ntype + 1
               IWK( lstype + iel ) = ntype
               IWK( lssiwt + ntype ) = liwfre
               IWK( lsiwtr + liwfre ) = ninvr
               IWK( lsiwtr + liwfre + 1 ) = nelvr
               IWK( lsswtr + ntype ) = lwfree
               liwfre = liwfre + 2 + ninvr + nelvr
               lwfree = lwfree + ninvr * nelvr
            ELSE
              IWK( lstype + iel ) = 0
            END IF
  350    CONTINUE

!  FOR EACH TYPE OF ELEMENT WITH INTERNAL VARIABLES:

         DO 360 i = 1, ntype
            liwfre = IWK( lssiwt + i )
            lwfree = IWK( lsswtr + i )
            ninvr = IWK( lsiwtr + liwfre )

!  FACTORIZE W. USE GAUSSIAN ELIMINATION WITH COMPLETE PIVOTING.
!  DETERMINE THE "MOST INDEPENDENT" SET OF COLUMNS OF W.

      CALL GELIM( ninvr, IWK( lsiwtr + liwfre + 1 ), &
                         IWK( lsiwtr + liwfre + 2 ), &
                         IWK( lsiwtr + liwfre + ninvr + 2 ), &
                         WK ( lswtra + lwfree ) )
  360    CONTINUE

!  COMPRESS THE DATA STRUCTURES TO REMOVE REDUNDANT INFORMATION.

         IF ( ntype < nel ) THEN
            k = lsswtr + ntype

!  COMPRESS INTEGER DATA.

            DO 370 i = 1, ntype
               IWK( k + i ) = IWK( lssiwt + i )
  370       CONTINUE
            lssiwt = k
         END IF
         k = lssiwt + ntype
         lniwtr = 0
         lnwtra = 0
         DO 400 i = 1, ntype
            liwfro = IWK( lssiwt + i ) - 1
            ninvr = IWK( lsiwtr + liwfro + 1 )
            DO 380 j = 1, 2 * ninvr + 2
               IWK( k + lniwtr + j ) = IWK( lsiwtr + liwfro + j )
  380       CONTINUE
            IWK( lssiwt + i ) = lniwtr + 1
            lniwtr = lniwtr + 2 + 2 * ninvr

!  COMPRESS REAL DATA.

            lwfreo = IWK( lsswtr + i ) - 1
            DO 390 j = 1, ninvr * ninvr
               WK( lswtra + lnwtra + j ) = WK( lswtra + lwfreo + j )
  390       CONTINUE
            IWK( lsswtr + i ) = lnwtra + 1
            lnwtra = lnwtra + ninvr * ninvr
  400    CONTINUE

!  RECORD THE LENGTHS OF THE PARTITIONS OF THE WORKSPACE USED.

         lsiwtr = k
         lwkstr = lswtra + lnwtra

! ---------------------------------------------------------------------
!  THE LIST OF VARIABLES is ALLOCATED TO nsets DISJOINTS SETS.
!  VARIABLE i OCCURS IN SET IWK( lsiset + i ).
! ---------------------------------------------------------------------

         lsiset = lsiwtr + lniwtr
         lssvse = lsiset + n
         lsend = lssvse + n + 1
         IF ( liwk < lsend + n ) THEN
            WRITE( iout, 2030 ) lsend + n - liwk
            inform = 4
            RETURN
         END IF

!  ASSIGN INITIAL SET NUMBERS TO EACH VARIABLE.

         nsets = 0
!DIR$    IVDEP
         DO 410 i = 1, n
            IWK( lsiset + i ) = n
  410    CONTINUE

!  USE THE CURTIS-POWELL-REID ALGORITHM TO DETERMINE WHICH SET EACH
!  VARIABLE BELONGS TO. LOOP OVER THE VARIABLES.

         DO 500 i = 1, n

!  LOOP OVER THE ELEMENTS WHICH USE VARIABLE I.
!  THE ELEMENTS ARE OBTAINED FROM A LINKED-LIST.

            VRUSED = .FALSE.
            ipt = IWK( lsptrs + i )
            IF ( ipt >= 0 ) THEN
               iell = IWK( lselts + i )
  420          CONTINUE
               iel = IELING( iell )
               itype = IWK( lstype + iel )
!              WRITE( 6, * ) ' ELEMENT ', iel

!  CHECK THAT THE VARIABLE BELONGS TO THE "INDEPENDENCE" SET OF
!  ELEMENTS WITH INTERNAL VARIABLES.

               IF ( itype > 0 ) THEN
                  liwfre = IWK( lssiwt + itype )
                  ninvr = IWK( lsiwtr + liwfre )
                  DO 430 j = 1, ninvr
                     k = j - 1
                     l = IWK( lsiwtr + liwfre + ninvr + 1 + j ) - 1
                     IF ( i == IELVAR( ISTAEV( iel ) + l ) ) &
                        GO TO 440
  430             CONTINUE
                  GO TO 450
  440             CONTINUE
               END IF
               VRUSED = .TRUE.
  450          CONTINUE

!  LOOP OVER THE COMPLETE LIST OF VARIABLES USED BY ELEMENT IEL.

!DIR$ IVDEP
               DO 460 j = ISTAEV( iel ), ISTAEV( iel + 1 ) - 1

!  IF VARIABLE IV is USED, FLAG THE SET THAT CONTAINS IT.

                  IWK( lsend + IWK( lsiset + IELVAR( j ) ) ) = 1
  460          CONTINUE

!  CHECK THE LINK-LIST TO SEE IF FURTHER ELEMENTS USE THE VARIABLE.

               IF ( ipt > 0 ) THEN
                  iell = IWK( lselts + ipt )
                  ipt = IWK( lsptrs + ipt )
                  GO TO 420
               END IF
            END IF

!  SEE IF THE VARIABLE MAY BE PLACED IN THE FIRST nsets SETS.

            IF ( VRUSED ) THEN
               DO 470 j = 1, nsets
                  IF ( IWK( lsend + j ) == 0 ) GO TO 480
                  IWK( lsend + j ) = 0
  470          CONTINUE

!  THE VARIABLE NEEDS A NEW SET.

               nsets = nsets + 1
               j = nsets

!  THE VARIABLE WILL BE PLACED IN SET J.

  480          CONTINUE
               IWK( lsiset + i ) = j

!  RESET THE FLAGS TO ZERO.

!DIR$ IVDEP
               DO 490 k = j, nsets
                  IWK( lsend + k ) = 0
  490          CONTINUE
            ELSE

!  THE VARIABLE is NOT TO BE USED.

               IWK( lsiset + i ) =  n
            END IF
  500    CONTINUE

!  CHECK THAT THERE is AT LEAST one SET.

         IF ( nsets /= 0 ) THEN

!  PRINT OUTPUT.

!DIR$ IVDEP
            DO 510 i = 1, n
               IWK( lsiset + i ) = MIN( IWK( lsiset + i ), nsets + 1 )
!              WRITE( 6, * ) ' VARIABLE ', i, ' SET ',
!    *                       IWK( lsiset + i )
  510       CONTINUE

! ---------------------------------------------------------------------
!  OBTAIN A LIST, IWK(LSEND), OF THE VARIABLES CORRESPONDING TO EACH SET.
! ---------------------------------------------------------------------

!  CLEAR IWK( lssvse ).

!DIR$ IVDEP
            DO 520 j = 2, nsets + 2
              IWK( lssvse + j ) = 0
  520       CONTINUE

!  COUNT THE NUMBER OF ELEMENTS IN EACH SET AND STORE IN IWK( lssvse ).
!  NEGATE THE SET NUMBERS IN IWK( lsiset ), SO THAT THEY ARE FLAGGED
!  AS IWK( lsiset ) is GRADUALLY OVERWRITTEN BY VARIABLE INDICES.

            DO 530 k = 1, n
               j = IWK( lsiset + k )
               IWK( lsiset + k ) = - j
               IWK( lssvse + j + 1 ) = IWK( lssvse + j + 1 ) + 1
  530       CONTINUE

!  COMPUTE THE STARTING ADDRESSES FOR EACH SET WITHIN IWK( lsiset ).

            IWK( lssvse + 1 ) = 1
            DO 540 j = 2, nsets + 2
               IWK( lssvse + j ) = IWK( lssvse + j ) + &
                                   IWK( lssvse + j - 1 )
  540       CONTINUE

!  STORE IN IWK( lsend ) THE VARIABLE WHOSE SET NUMBER
!  is THE IWK( lssvse + j )-TH ENTRY OF IWK( lsend ).

            isofar = 0
            DO 570 j = 1, nsets + 1
               istarj = IWK( lssvse + j )
               DO 550 ivarp1 = isofar + 1, n
                  IF ( istarj < ivarp1 ) GO TO 560
  550          CONTINUE
               ivarp1 = n + 1
  560          CONTINUE
               isofar = ivarp1 - 1
               IWK( lsend + j ) = isofar
  570       CONTINUE

!  REORDER THE ELEMENTS INTO SET ORDER.
!  FILL IN EACH SET FROM THE FRONT. AS A NEW ENTRY is PLACED
!  IN SET k INCREASE THE POINTER IWK( lssvse + k ) BY one AND FIND
!  THE NEW VARIABLE, IWK( lsend + k ), THAT CORRESPONDS TO THE SET NOW
!  POINTED TO BY IWK( lssvse + k ).

            DO 660 j = 1, nsets + 1

!  DETERMINE THE next UNPLACED ENTRY, IWK( lssvse ), IN IWK( lsiset ).

  610          CONTINUE
               istrt = IWK( lssvse + j )

!  SEE IF ALL THE ELEMENTS IN SET j HAVE BEEN ASSIGNED.

               IF ( istrt == IWK( lssvse + j + 1 ) ) GO TO 660
               IF ( IWK( lsiset + istrt ) > 0 ) GO TO 660

!  EXTRACT THE VARIABLE AND SET NUMBERS OF THE STARTING ELEMENT.

               ivar = IWK( lsend + j )
               jset = - IWK( lsiset + istrt )

!  MOVE ELEMENTS IN A CYCLE, ENDING BACK AT SET J.

               DO 640 k = istrt, n

!  FIND THE FIRST EMPTY LOCATION IN SET jset IN IWK( lsend )

                  inext = IWK( lssvse + jset )

!  EXTRACT THE VARIABLE INDEX OF THE next ELEMENT.

                  newvar = IWK( lsend + jset )

!  UPDATE IWK( lssvse + jset ), FIND THE NEW VARIABLE INDEX AND STORE
!  IT IN IWK( lsend + jset ).

                  istarj = inext + 1
                  IWK( lssvse + jset ) = istarj
                  DO 620 ivarp1 = newvar + 1, n
                     IF ( istarj < ivarp1 ) GO TO 630
  620             CONTINUE
                  ivarp1 = n + 1
  630             CONTINUE
                  IWK( lsend + jset ) = ivarp1 - 1

!  IF THE ENTRY BELONGS IN THE J-TH SET, THE CYCLE is COMPLETE.

                  IF ( jset == j ) GO TO 650

!  EXTRACT THE NUMBER OF THE SET OF THE next ELEMENT.

                  newset = - IWK( lsiset + inext )

!  STORE THE VARIABLE INDEX OF THE CURRENT ELEMENT.

                  IWK( lsiset + inext ) = ivar

!  MAKE THE next ELEMENT INTO THE CURRENT ONE.

                  ivar = newvar
                  jset = newset
  640          CONTINUE

!  THE CYCLE is COMPLETE.

  650          CONTINUE

!  STORE THE VARIABLE INDEX OF THE STARTING ELEMENT.

               IWK( lsiset + istrt ) = ivar
               GO TO 610
  660       CONTINUE

!  REVISE IWK( lssvse ) TO POINT TO THE START OF EACH SET.

            DO 670 j = nsets + 1, 2, - 1
               IWK( lssvse + j ) = IWK( lssvse + j - 1 )
  670       CONTINUE
            IWK( lssvse + 1 ) = 1
         END IF
!        DO 671 i = 1, nsets
!           WRITE( 6, * ) ' SET ', i, ' VARIABLES ',
!    * ( IWK( lsiset + j ), j = IWK( lssvse + i ),
!    *         IWK( lssvse + i + 1 ) - 1 )
! 671    CONTINUE
      ELSE

!  EXACT GRADIENTS ARE USED. NO FURTHER PARTITIONING OF THE WORKSPACE
!  is NEEDED.

         lstype = lsend
         lsswtr = lstype
         lssiwt = lsswtr
         lsiwtr = lssiwt
         lsiset = lsiwtr
         lssvse = lsiset
         lswtra = lwkstr
      END IF

!  SET THE LENGTH OF THE REMAINING PARTITIONS OF THE WORKSPACE FOR
!  ARRAY BOUND CHECKING IN CALLS TO OTHER SUBPROGRAMS.

      lntype = MAX( 1, lsswtr - lstype )
      lnswtr = MAX( 1, lssiwt - lsswtr )
      lnsiwt = MAX( 1, lsiwtr - lssiwt )
      lniwtr = MAX( 1, lsiset - lsiwtr )
      lniset = MAX( 1, lssvse - lsiset )
      lnsvse = MAX( 1, lsend - lssvse )
      lnwtra = MAX( 1, lwkstr - lswtra )

!  RECORD THE LENGTHS OF THE REMAINING INTEGER AND REAL WORKSPACE.

      liwk2 = liwk - lsend
      lwk2 = lwk - lwkstr

! -- SET THE STARTING ADDRESSES FOR THE PARTITIONS WITHIN FUVALS. --

!  A FULL DESCRIPTION OF THE PARTITIONS OF FUVALS is GIVEN IN THE
!  THE INTRODUCTORY COMMENTS TO SUBROUTINE SBMIN.

      lfxi = 0
      lgxi = lfxi + nel
      lhxi = INTVAR( nel1 ) - 1
      lggfx = lggfx - 1
      ldx = lggfx + n
      lgrjac = ldx + n
      lend = lgrjac + nvargp

!  PRINT ALL OF THE STARTING ADDRESSES FOR THE WORKSPACE ARRAY
!  PARTITIONS.

      IF ( iprint >= 3 ) THEN
         WRITE( iout, 2000 ) &
           lfxi, lgxi, lhxi, lggfx, ldx, lgrjac, lend, lfuval
         WRITE( iout, 2010 ) lsptrs, lselts, lindex, &
                             lswksp, liused, lfreec, lnnonz, lnonz2, &
                             lsymmd, lsymmh, lslgrp, lstajc, &
                             lstagv, lsvgrp, lgcolj, lvaljr, &
                             lstype, lsswtr, lssiwt, lsiwtr, lsiset, &
                             lssvse, lsend, liwk
         WRITE( iout, 2020 ) lqgrad, lbreak, lp,     lxcp, lx0, lgx0, &
                             ldeltx, lbnd, lswtra, lwkstr, lwk
      END IF

!  CHECK THAT THE ARRAY FUVALS HAS SUFFICIENT ROOM FOR THE CALCULATION.

      IF ( lend > lfuval ) THEN
         WRITE( iout, 2050 ) lend - lfuval
         inform = 5
         RETURN
      END IF

!  SET THE LENGTH OF EACH PARTITION OF THE REAL WORKSPACE ARRAY
!  FUVALS FOR ARRAY BOUND CHECKING IN CALLS TO OTHER SUBPROGRAMS.

      lnfxi = MAX( 1, lgxi - lfxi )
      lngxi = MAX( 1, lhxi - lgxi )
      lnguvl = MAX( 1, lhxi - lfxi )
      lnhxi = MAX( 1, lggfx - lhxi )
      lnhuvl = MAX( 1, lggfx - lfxi )
      lnggfx = MAX( 1, ldx - lggfx )
      lndx = MAX( 1, lgrjac - ldx )
      lngrjc = MAX( 1, lend - lgrjac )
      maxsin = MAX( 1, maxsin )
      maxsel = MAX( 1, maxsel )
      inform = 0

!  NON-EXECUTABLE STATEMENTS.

 2000 FORMAT( /, ' Starting addresses for the partitions of FUVALS ', &
              /, ' ----------------------------------------------- ', &
             //, '   lfxi   lgxi   lhxi  lggfx ', &
                 '   ldx lgrjac   lend   lfuval ', /, 7I7, I9 )
 2010 FORMAT( /, ' Starting addresses for partitions of IWK ', &
              /, ' ---------------------------------------- ', //, &
        ' lsptrs lselts lindex lswksp liused lfreec lnnonz lnonz2  ...', &
         /, 8I7, //, &
         ' ...... lsymmd lsymmh lslgrp lstajc lstagv lsvgrp lgcolj ...', &
         /, 7X, 7I7, //, &
        ' ...... lvaljr lstype lsswtr lssiwt lsiwtr lsiset lssvse ... ', &
         /, 7X, 7I7, //, &
         ' ......  lsend     liwk ', &
         /, 7X, I7, I9 )
 2020 FORMAT( /, ' Starting addresses for partitions of WK ', &
              /, ' --------------------------------------- ', //, &
         '   lqgrad   lbreak       lp     lxcp      lx0     lgx0 ... ', &
         /, 6I9, //, &
         '   ......   ldeltx     lbnd   lswtra   lwkstr      lwk ', &
         /, 9X, 5I9 )
 2030 FORMAT( /, ' INITW: The size of array IWK must be increased', &
                 ' by at least ', I12 )
 2040 FORMAT( /, ' INITW: The size of array WK must be increased', &
                 ' by at least ', I8 )
 2050 FORMAT( /, ' INITW: The size of array FUVALS must be increased', &
                 ' by at least ', I8 )
      RETURN

!  END OF SUBROUTINE INITW.

      END
