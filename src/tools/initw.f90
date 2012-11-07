! THIS VERSION: CUTEST 1.0 - 05/11/2012 AT 13:45 GMT.

!-*-*-*-*-*-*-  C U T E S T    C S E T U P    S U B R O U T I N E  -*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 77 version originally released in CUTE, 20th june 1990
!   fortran 2003 version released in CUTEst, 5th November 2012

      SUBROUTINE INITW( n, ng, nel, IELING, leling, ISTADG, lstadg,            &
                        IELVAR, lelvar, ISTAEV, lstaev, INTVAR, lntvar,        &
                        ISTADH, lstadh, ICNA, licna, ISTADA, lstada,           &
                        ITYPEE, litype, gxeqx, lgxeqx, intrep, lintre,         &
                        lfuval, ALTRIV, direct, FDGRAD, lfxi, lgxi,            &
                        lhxi, lggfx, ldx, lgrjac, lqgrad, lbreak, lp,          &
                        lxcp, lx0, lgx0, ldeltx, lbnd, lwkstr, lsptrs,         &
                        lselts, lindex, lswksp, lstagv, lstajc, liused,        &
                        lfreec, lnnonz, lnonz2, lsymmd, lsymmh, lslgrp,        &
                        lsvgrp, lgcolj, lvaljr, lsend, lnptrs, lnelts,         &
                        lnndex, lnwksp, lnstgv, lnstjc, lniuse, lnfrec,        &
                        lnnnon, lnnno2, lnsymd, lnsymh, lnlgrp, lnvgrp,        &
                        lngclj, lnvljr, lnqgrd, lnbrak, lnp, lnbnd,            &
                        lnfxi, lngxi, lnguvl, lnhxi, lnhuvl, lnggfx,           &
                        lndx, lngrjc, liwk2, lwk2, maxsin, ninvar,             &
                        ntype, nsets, maxsel, lstype, lsswtr, lssiwt,          &
                        lsiwtr, lswtra, lntype, lnswtr, lnsiwt, lniwtr,        &
                        lnwtra, lsiset, lssvse, lniset, lnsvse, RANGE,         &
                        IWK, liwk, WK, lwk, iprint, iout, status )
                        INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

! ------------------------------------------------------------------------
!  compute the starting addresses for the partitions of the workspace
!  arrays FUVALS, IWK AND WK. Also fill relevant portions of IWK. The
!  addresses for IWK and WK are described here; those for FUVALS are
!  described in the introductory comments to the LANCELOT subroutine SBMIN
! ------------------------------------------------------------------------

      INTEGER :: n, ng, nel, status, litype, lwk, liwk
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
      LOGICAL :: direct, ALTRIV, FDGRAD
      INTEGER :: IELVAR( lelvar ), ISTAEV( lstaev )
      INTEGER :: ISTADH( lstadh ), IWK ( liwk )
      INTEGER :: INTVAR( lntvar ), ISTADG( lstadg )
      INTEGER :: ICNA ( licna ), ISTADA( lstada )
      INTEGER :: IELING( leling ), ITYPEE( litype )
      REAL ( KIND = wp ) :: WK ( lwk )
      LOGICAL :: gxeqx ( lgxeqx ), INTREP( lintre )

!  Local variables

      INTEGER :: i, j, k, l, nvargp, iielts, ientry, ig, is
      INTEGER :: lend, nsizeh, nel1, ntotel, ng1, iel
      INTEGER :: lwfree, liwfre, nelvr, lw1, liwfro, lwfreo
      INTEGER :: iell, itype, isofar, istarj, ivarp1, IVAR
      INTEGER :: jset, inext, newvar, newset, ipt, istrt
      INTEGER :: ninvr, ii, jj, kk, ll
      LOGICAL :: NONTRV, ALLLIN, VRUSED
      REAL ( KIND = wp ) :: zero, one
      PARAMETER ( zero = 0.0_wp, one = 1.0_wp )

!  external subroutines and functions used

      EXTERNAL :: SYMMH

!  set constants

      nel1 = nel + 1
      ng1 = ng + 1
      ntotel = ISTADG( ng1 ) - 1
      ALLLIN = nel == 0

!  set up INTVAR, the starting addresses for the element gradients with respect
!  to their internal variables. Also compute maxsin, the maximum number of 
!  internal variables in an element

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

!  compute the total number of internal variables

      ninvar = INTVAR( nel1 ) - INTVAR( 1 )

!  calculate the length, iielts, of workspace required to determine which 
!  elements use each of the variables. Also find the maximum number of 
!  variables in an element, maxsel. This is a dummy run for loop 130 merely 
!  to calculate the space required

      IF ( liwk < n ) THEN
         WRITE( iout, 2030 ) n - liwk
         status = 4
         RETURN
      END IF

!  IWK(i) will be used as a list of links chaining the elements using
!  variable i. If IWK(i) is negative, the list is empty

!DIR$ IVDEP
      DO 20 i = 1, n
         IWK( i ) = - 1
   20 CONTINUE
      iielts = n
      maxsel = 0
      IF ( .NOT. ALLLIN ) THEN

!  loop over the group, considering each nonlinear element in turn

         DO 50 i = 1, ntotel
            iel = IELING( i )
            maxsel = MAX( maxsel, ISTAEV( iel + 1 ) - ISTAEV( iel ) )

!  loop on the variables from the i-th element

            DO 40 k = ISTAEV( iel ), ISTAEV( iel + 1 ) - 1
               ientry = IELVAR( k )
               IF ( IWK( ientry ) >= 0 ) THEN

!  if we have reached the end of the list of the elements using the variable 
!  IELVAR(k), add the iel-th element to it. otherwise, find the next entry 
!  in the list

   30             CONTINUE
                  IF ( IWK( ientry ) > 0 ) THEN
                     ientry = IWK( ientry )
                     GO TO 30
                  ELSE
                     iielts = iielts + 1
                     IF ( iielts > liwk ) THEN
                        WRITE( iout, 2030 ) iielts - liwk
                        status = 4
                        RETURN
                     END IF
                     IWK( ientry ) = iielts
                     IWK( iielts ) = 0
                  END IF
               ELSE

!  the list of elements involving the variable ielvar( k ) was previously 
!  empty. indicate that the list has now been started and that its end has 
!  been reached

                  IWK( ientry ) = 0
               END IF
   40       CONTINUE
   50    CONTINUE
      END IF

! -- calculate the starting addresses for the integer workspace --

!  IWK(lsptrs + j), j = 1, ..., iielts,  will contain the links for the lists 
!  of nonlinear elements which use each variable

      lsptrs = 0

!  IWK(lselts + j), j = 1, ..., iielts, will contain the lists of nonlinear 
!  elements corresponding to the previously mentioned links

      lselts = lsptrs + iielts

!  IWK(lindex + j), j = 1, ..., n, will contain the status of the
!  j-th variable as the current iteration progresses. Possible values
!  are 0 if the variable lies away from its bounds, 1 and 2 if it lies
!  on its lower or upper bounds (respectively) - these may be problem
!  bounds or trust region bounds, and 3 if the variable is fixed.

      lindex = lselts + iielts

!  IWK( lswksp + j ), j = 1, ..., MAX( ntotel, nel, n + n ), is used for
!  workspace by the matrix-vector product subroutine HESPRD

      lswksp = lindex + n

!  IWK( liused + j ), j = 1, ..., MAX( n, ng ) will be used as workspace by 
!  the matrix-vector product subroutine HESPRD

      IF ( direct ) THEN
         liused = lswksp + MAX( ntotel, nel, n + n )
      ELSE
         liused = lswksp + MAX( ntotel, nel, n )
      END IF

!  IWK(lfreed + j), j = 1, ..., nfreec will give the indices of the variables 
!  that are considered to be free from their bounds at the current generalized 
!  Cauchy point

      lfreec = liused + MAX( n, ng )

!  IWK(lnnonz + j), j = 1, ..., nnnonz will give the indices of the nonzeros 
!  in the vector obtained as a result of the matrix-vector product from 
!  subroutine HESPRD

      lnnonz = lfreec + n

!  IWK(lnonz2 + j), j = 1, ..., ng, will be used as further workspace by the 
!  matrix-vector product subroutine HESPRD

      lnonz2 = lnnonz + n

!  IWK(lsymmd + j), j = 1, ..., maxsin, will give the location of the j-th 
!  diagonal of a maxsin by maxsin symmetric matrix in an upper triangular 
!  storage scheme

      lsymmd = lnonz2 + ng

!  IWK(lsymmh + ( i - 1 ) * maxsin + j), i, j = 1, ..., maxsin, will give the 
!  location of the (i,j)-th entry of a maxsin by maxsin symmetric matrix in 
!  an upper triangular storage scheme.

      lsymmh = lsymmd + maxsin

!  IWK( lslgrp + j), j = 1, ..., ntotel, will contain the number of the group
!  that uses nonlinear element j

      lslgrp = lsymmh + maxsin * maxsin

!  IWK(lstajc + j), j = 1, ..., n, will contain the starting addresses for the 
!  list of nontrivial group  which use the j-th variable. iwk(lstajc + n + 1) 
!  will point to the first free location in iwk after the list of nontrivial 
!  group  for the n-th variable

      lstajc = lslgrp + ntotel

!  IWK(lstagv + j), j = 1, ..., ng, will contain the starting addresses for 
!  the list of variables which occur in the j-th group. IWK(lstagv + ng + 1) 
!  will point to the first free location in IWK after the list of variables 
!  for the ng-th group

      lstagv = lstajc + n + 1

!  IWK(lsvgrp + j), j = 1, ..., nvargp, will contain the indices of the 
!  variables which are used by each group in turn. those for group i occur in 
!  locations IWK(lstagv + i) to iwk(lstagv + i + 1) - 1

      lsvgrp = lstagv + ng + 1

!  check that there is sufficient workspace.

      IF ( lsvgrp > liwk ) THEN
         WRITE( iout, 2030 ) lsvgrp - liwk
         status = 4
         RETURN
      END IF

!  IWK(lgcolj + j), j = 1, ..., nvargp, will contain the indices of the
!  nontrivial group  which use each variable in turn. those for variable
!  i occur in locations IWK(lstajc + i) TO  IWK(lstajc + i + 1) - 1

      lgcolj = lsvgrp + 1

!  determine which elements use each variable. Initialization

      IF ( .NOT. alllin ) THEN

!  IWK(i) will be used as a list of links chaining the elements using
!  variable i. If IWK(i) is negative, the list is empty

!DIR$ IVDEP
         DO 100 i = 1, n
            IWK( i ) = - 1
  100    CONTINUE
         iielts = n

!  loop over the group, considering each nonlinear element in turn

         DO 130 i = 1, ntotel
            iel = IELING( i )

!  loop on the variables of the i-th element

            DO 120 k = ISTAEV( iel ), ISTAEV( iel + 1 ) - 1
               ientry = IELVAR( k )
               IF ( IWK( ientry ) >= 0 ) THEN

!  If we have reached the end of the list of the elements using the variable 
!  IELVAR(k), add the i-th element to it and record that the end of the list 
!  has occured. Otherwise, find the next entry in the list

  110             CONTINUE
                  IF ( IWK( ientry ) > 0 ) THEN
                     ientry = IWK( ientry )
                     GO TO 110
                  ELSE
                     iielts = iielts + 1
                     IWK( ientry ) = iielts
                     IWK( lselts + iielts ) = i
                     IWK( iielts ) = 0
                  END IF
               ELSE

!  the list of elements involving the variable IELVAR(k) was previously empty.
!  Indicate that the list has now been started, record the element that 
!  contains the variable and indicate that the end of the list has been reached

                  IWK( lselts + ientry ) = i
                  IWK( ientry ) = 0
               END IF
  120       CONTINUE
  130    CONTINUE
      END IF

!  set up symmetric addresses for the upper triangular storage schemes for the 
!  element Hessians

      IF ( maxsin > 0 ) CALL SYMMH( maxsin, IWK( lsymmh + 1 ),                 &
                                    IWK( lsymmd + 1 ) )

!  set up the starting addresses for the element Hessians with respect to 
!  their internal variables and a pointer beyond the end of the space required 
!  for the Hessians

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

!  WK(lqgrad + j), j = 1, ..., n, WILL CONTAIN THE GRADIENT OF
!  THE QUADRATIC MODEL AT THE CURRENT ESTIMATE OF THE MINIMIZER

      lqgrad = MAX( MAX( ng, maxsel ) + 2 * maxsin, &
                    n + ninvar + MAX( maxsel, ninvar ) )

!  WK(lbreak + j), j = 1, ..., n, WILL CONTAIN THE BREAKPOINTS
!  ALONG THE CAUCHY ARC FROM THE CURRENT ESTIMATE OF THE MINIMIZER

      lbreak = lqgrad + n

!  WK(lp + j), j = 1, ..., n, WILL CONTAIN THE VECTOR REQUIRED
!  BY THE MATRIX-VECTOR PRODUCT SUBROUTINE HESPRD

      lp = lbreak + n

!  WK(lxcp + j), j = 1, ..., n, WILL CONTAIN THE CURRENT
! (APPROXIMATE) CAUCHY POINT

      lxcp = lp + n

!  WK(lx0 + j), j = 1, ..., n, WILL CONTAIN THE CURRENT
!  START OF THE CAUCHY SEARCH. THIS FEATURE is ONLY USED IF
!  MORE THAN one CYCLE is USED TO SOLVE THE BQP ACCURATELY

      lx0 = lxcp + n

!  WK(lgx0 + j), j = 1, ..., n, WILL CONTAIN THE GRADIENT AT THE
!  CURRENT START OF THE CAUCHY SEARCH. THIS FEATURE is ONLY USED IF
!  MORE THAN one CYCLE is USED TO SOLVE THE BQP ACCURATELY

      lgx0 = lx0 + n

!  WK(ldeltx + j), j = 1, ..., n, WILL CONTAIN THE STEP TAKEN
!  DURING A SINGLE CYCLE, IF MORE THAN one CYCLE is USED TO SOLVE
!  THE BQP ACCURATELY

      ldeltx = lgx0 + n

!  WK(lbnd + 2 * ( j - 1 ) + i), j = 1, ..., n, i = 1, 2, WILL
!  CONTAIN THE CURRENT LOWER (I=1) AND UPPER (I=2) BOUNDS ON THE
!  VARIABLES DEFINED BY THE INTERSECTION OF THE TRUST REGION WITH
!  THE FEASIBLE BOX

      lbnd = ldeltx + n

!  WK(lwkstr + j), j = 1, ..., lwk2, is THE REMAINING REAL
!  WORKSPACE WHICH is FREE FOR OTHER PURPOSES, SUCH AS FOR FORMING
!  THE FACTORIZATION OF THE MODEL HESSIAN, IF REQUIRED

      lwkstr = lbnd + n + n
      lwk2 = lwk - lwkstr

!  check that there is sufficient real workspace

      IF ( lwkstr > lwk ) THEN
         status = 5
         WRITE( iout, 2040 ) lwkstr - lwk
         RETURN
      END IF

!  set the length of each partition of the real workspace array
!  fuvals for array bound checking in calls to other subprograms

      lnqgrd = lbreak - lqgrad
      lnbrak = lp - lbreak
      lnp = lxcp - lp
      lnbnd = lwkstr - lbnd

!  store the indices of variables which appears in each group
!  and how many group  use each variable. start by initializing
!  counting arrays to zero

!dir$ ivdep
      DO 150 j = 1, n
         IWK( lswksp + j ) = 0
         IWK( lstajc + j + 1 ) = 0
  150 CONTINUE

!  altriv specifies whether all the group  are trivial

      ALTRIV = .TRUE.

!  count the total number of variables in all the group, nvargp

      nvargp = 0
      IWK( lstagv + 1 ) = 1

!  loop over the group . see if the ig-th group is trivial

      DO 200 ig = 1, ng
         NONTRV = .NOT. GXEQX( ig )

!  check to see if all of the group  are trivial

         IF ( NONTRV ) ALTRIV = .FALSE.

!  loop over the nonlinear elements from the ig-th group

         DO 170 k = ISTADG( ig ), ISTADG( ig + 1 ) - 1
            iel = IELING( k )

!  run through all the elemental variables changing the i-th entry of
!  iwk( lswksp ) from zero to one if variable i appears in an element

            DO 160 j = ISTAEV( iel ), ISTAEV( iel + 1 ) - 1
               i = IELVAR( j )
               IF ( IWK( lswksp + i ) == 0 ) THEN
                  IWK( lswksp + i ) = 1

!  if there is sufficient room, record the nonlinear variables from
!  the ig-th group

                  IF ( lgcolj > liwk ) THEN
                     WRITE( iout, 2030 ) lgcolj - liwk
                     status = 4
                     RETURN
                  END IF
                  IWK( lgcolj ) = i
                  lgcolj = lgcolj + 1
                  nvargp = nvargp + 1
               END IF
  160       CONTINUE

!  record that nonlinear element k occurs in group ielgrp(iel)

            IWK( lslgrp + k ) = ig
  170    CONTINUE

!  consider variables which arise from the linear element

         DO 180 j = ISTADA( ig ), ISTADA( ig + 1 ) - 1
            i = ICNA( j )
            IF ( IWK( lswksp + i ) == 0 ) THEN
               IWK( lswksp + i ) = 1

!  if there is sufficient room, record the linear variables from
!  the ig-th group

               IF ( lgcolj > liwk ) THEN
                  WRITE( iout, 2030 ) lgcolj - liwk
                  status = 4
                  RETURN
               END IF
               IWK( lgcolj ) = i
               lgcolj = lgcolj + 1
               nvargp = nvargp + 1
            END IF
  180    CONTINUE

!  reset the status array iwk(lswksp) to zero

         DO 190 j = IWK( lstagv + ig ), nvargp
            l = IWK( lsvgrp + j )
            IWK( lswksp + l ) = 0

!  record that one further nontrivial group uses variable l

            IF ( NONTRV ) IWK( lstajc + l + 1 ) = &
                          IWK( lstajc + l + 1 ) + 1
  190    CONTINUE

!  record the starting address of the variables in the next group.

         IWK( lstagv + ig + 1 ) = nvargp + 1
  200 CONTINUE

! -- continue setting starting addresses for partitions of iwk --

!  iwk(lvaljr + j), j = 1, ..., nvargp, will contain the positions in
!  fuvals (relative to the starting address lgrjac) of the nonzeros of
!  the jacobian of the group  corresponding to the variables as ordered
!  in iwk( lsvgrp + j ).

      lvaljr = lgcolj + nvargp

!  lsend gives the total fixed amount of integer workspace used. liwk2
!  gives the amount of workspace which is free for other purposes, such
!  as for forming the factorization of the model hessian, if required

      lsend = lvaljr + nvargp

!  check that there is sufficient integer workspace.

      IF ( lsend > liwk ) THEN
         WRITE( iout, 2030 ) lsend - liwk
         status = 4
         RETURN
      END IF

!  set the length of each partition of the integer workspace for
!  array bound checking in calls to other subprograms.

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

!  set the starting addresses for the lists of nontrivial group 
!  which use each variable in turn

      k = lstajc + 1
      IWK( k ) = 1
      DO 210 i = 2, n + 1
         k = k + 1
         IWK( k ) = IWK( k ) + IWK( k - 1 )
  210 CONTINUE

!  consider the ig-th group in order to associate variables with group

      DO 230 ig = 1, ng
         IF ( .NOT. GXEQX( ig ) ) THEN
            DO 220 i = IWK( lstagv + ig ), IWK( lstagv + ig + 1 ) - 1
                l = lstajc + IWK( lsvgrp + i )

!  record that group ig uses variable iwk(lsvgrp + i)

               j = IWK( l )
               IWK( lgcolj + j ) = ig

!  store the locations in the jacobian of the group  of the nonzeros
!  corresponding to each variable in the ig-th group. increment the
!  starting address for the pointer to the next group using variable
!  iwk(lsvgrp + i)

               IWK( lvaljr + i ) = j
               IWK( l ) = j + 1
  220      CONTINUE
         END IF
  230 CONTINUE

!  reset the starting addresses for the lists of group  using
!  each variable.

      DO 240 i = n, 2, - 1
         l = lstajc + i
         IWK( l ) = IWK( l - 1 )
  240 CONTINUE
      IWK( lstajc + 1 ) = 1

!  initialize workspace values for subroutine hesprd

!dir$ ivdep
      DO 250 j = 1, MAX( n, ng )
         IWK( liused + j ) = 0
  250 CONTINUE

!  define further partitions of the workspace whenever finite-
! -difference gradients are used

      IF (  FDGRAD ) THEN

! -- continue setting starting addresses for partitions of iwk --

!  the range transformation for each nonlinear element is of a given
!  type. suppose there are ntype non-trivial types. iwk( lstype + i )
!  gives the type of nonlinear element i for i = 1, ...., nel

         lstype = lsend

!  the range transformation from elemental to internal variables is
!  defined by a matrix w. for each non-trivial transformation, the
!  matrix w is recorded. the information for the i-th type starts in
!  location lswtra + iwk( lsswtr + i ), i = 1, ...., ntype

         lsswtr = lstype + nel

!  for each type of nonlinear element using a nontrivial range
!  transformation, integer information isd also recorded.
!  the information for the i-th type starts in location
!  lsiwtr + iwk( lssiwt + i ), i = 1, ...., ntype

         lssiwt = lsswtr + nel

!  the following pieces of integer information are recorded about
!  the i-th type of nonlinear element:

!    iwk( lsiwtr + iwk( lssiwt + i ) + 1 ):
!            the number of internal variables, ninvr.
!    iwk( lsiwtr + iwk( lssiwt + i ) + 2 ):
!            the number of elemental variables, nelvr.
!    iwk( lsiwtr + iwk( lssiwt + i ) + 2 + i ),
!         i = 1, ..., nelvr + ninvr:
!            pivot sequences for the lu factors of w

!  after the factorization and compression, only ninvr linearly
!  independent columns of w are stored

         lsiwtr = lssiwt + nel
         IF ( liwk < lsiwtr ) THEN
            WRITE( iout, 2030 ) lsiwtr - liwk
            status = 4
            RETURN
         END IF

! -- continue setting starting addresses for partitions of wk --

!  the following pieces of integer information are recorded about
!  the i-th type of nonlinear element:

!    iwk( lswtra + iwk( lsswtr + i ) + i ),
!         i = 1, ..., nelvr * ninvr:  the matrix w stored by columns

!  after the factorization and compression, only ninvr linearly
!  independent columns of w are stored

         lswtra = lwkstr

! ---------------------------------------------------------------------
!  consider only elements which use internal variables
! ---------------------------------------------------------------------

         ntype = 0
         lwfree = 1
         liwfre = 1

!  loop over all nonlinear elements

         DO 350 iel = 1, nel
            IF ( INTREP( iel ) ) THEN

!  calculate the range transformation matrix w

               is = ISTAEV( iel )
               ninvr = INTVAR( iel + 1 ) - INTVAR( iel )
               nelvr = ISTAEV( iel + 1 ) - is
               lw1 = lwfree + ninvr * nelvr - 1
               l = lswtra + lw1
!dir$ ivdep
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

!  check to see if any of the columns belong to duplicated variables

                  ii = IELVAR( is + i )
                  DO 290 j = 1, i - 1
                     IF ( IELVAR( is + j ) == ii ) GO TO 300
  290             CONTINUE
                  GO TO 320

!  amalgamate columns from duplicate variables

  300             CONTINUE
                  kk = lswtra + lwfree + ( j - 1 ) * ninvr - 1
                  ll = k - ninvr - 1
                  DO 310 jj = 1, ninvr
                     WK( kk + jj ) = WK( kk + jj ) + WK( ll + jj )
                     WK( ll + jj ) = zero
  310             CONTINUE
  320          CONTINUE

!  compare this transformation matrix with previous ones

               DO 340 i = 1, ntype
                  IF ( IWK( lsiwtr + IWK( lssiwt + i ) ) /= ninvr &
                       .OR. IWK( lsiwtr + IWK( lssiwt + i ) + 1 ) &
                       /= nelvr ) GO TO 340
                  DO 330 j = 0, ninvr * nelvr - 1
                     IF ( WK( lswtra + lwfree + j ) /= WK( &
                          lswtra + IWK( lsswtr + i ) + j ) ) GO TO 340
  330             CONTINUE

!  the transformation is an existing one. record which one

                  IWK( lstype + iel ) = i
                  GO TO 350
  340          CONTINUE

!  ensure that there is sufficient room

               IF ( liwfre + 2 + ninvr + nelvr > liwk ) THEN
                  status = 4
                  WRITE( iout, 2030 ) liwfre + 2 + ninvr + nelvr - liwk
                  RETURN
               END IF
               IF ( lwfree + ninvr * nelvr > lwk ) THEN
                  status = 5
                  WRITE( iout, 2040 ) lwfree + ninvr * nelvr - lwk
                  RETURN
               END IF

!  the transformation defines a new type. record its details

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

!  for each type of element with internal variables:

         DO 360 i = 1, ntype
            liwfre = IWK( lssiwt + i )
            lwfree = IWK( lsswtr + i )
            ninvr = IWK( lsiwtr + liwfre )

!  factorize w. use gaussian elimination with complete pivoting.
!  determine the "most independent" set of columns of w

      CALL GELIM( ninvr, IWK( lsiwtr + liwfre + 1 ), &
                         IWK( lsiwtr + liwfre + 2 ), &
                         IWK( lsiwtr + liwfre + ninvr + 2 ), &
                         WK ( lswtra + lwfree ) )
  360    CONTINUE

!  compress the data structures to remove redundant information

         IF ( ntype < nel ) THEN
            k = lsswtr + ntype

!  compress integer data

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

!  compress real data

            lwfreo = IWK( lsswtr + i ) - 1
            DO 390 j = 1, ninvr * ninvr
               WK( lswtra + lnwtra + j ) = WK( lswtra + lwfreo + j )
  390       CONTINUE
            IWK( lsswtr + i ) = lnwtra + 1
            lnwtra = lnwtra + ninvr * ninvr
  400    CONTINUE

!  record the lengths of the partitions of the workspace used

         lsiwtr = k
         lwkstr = lswtra + lnwtra

! ---------------------------------------------------------------------
!  the list of variables is allocated to nsets disjoints sets.
!  variable i occurs in set iwk( lsiset + i )
! ---------------------------------------------------------------------

         lsiset = lsiwtr + lniwtr
         lssvse = lsiset + n
         lsend = lssvse + n + 1
         IF ( liwk < lsend + n ) THEN
            WRITE( iout, 2030 ) lsend + n - liwk
            status = 4
            RETURN
         END IF

!  assign initial set numbers to each variable

         nsets = 0
!dir$    ivdep
         DO 410 i = 1, n
            IWK( lsiset + i ) = n
  410    CONTINUE

!  use the curtis-powell-reid algorithm to determine which set each
!  variable belongs to. loop over the variables

         DO 500 i = 1, n

!  loop over the elements which use variable i. the elements are obtained from 
!  a linked-list

            vrused = .FALSE.
            ipt = IWK( lsptrs + i )
            IF ( ipt >= 0 ) THEN
               iell = IWK( lselts + i )
  420          CONTINUE
               iel = IELING( iell )
               itype = IWK( lstype + iel )
!              write( 6, * ) ' element ', iel

!  check that the variable belongs to the "independence" set of
!  elements with internal variables

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
               vrused = .TRUE.
  450          CONTINUE

!  loop over the complete list of variables used by element iel

!dir$ ivdep
               DO 460 j = ISTAEV( iel ), ISTAEV( iel + 1 ) - 1

!  if variable iv is used, flag the set that contains it

                  IWK( lsend + IWK( lsiset + IELVAR( j ) ) ) = 1
  460          CONTINUE

!  check the link-list to see if further elements use the variable

               IF ( ipt > 0 ) THEN
                  iell = IWK( lselts + ipt )
                  ipt = IWK( lsptrs + ipt )
                  GO TO 420
               END IF
            END IF

!  see if the variable may be placed in the first nsets sets

            IF ( VRUSED ) THEN
               DO 470 j = 1, nsets
                  IF ( IWK( lsend + j ) == 0 ) GO TO 480
                  IWK( lsend + j ) = 0
  470          CONTINUE

!  the variable needs a new set

               nsets = nsets + 1
               j = nsets

!  the variable will be placed in set j

  480          CONTINUE
               IWK( lsiset + i ) = j

!  reset the flags to zero

!dir$ ivdep
               DO 490 k = j, nsets
                  IWK( lsend + k ) = 0
  490          CONTINUE
            ELSE

!  the variable is not to be used

               IWK( lsiset + i ) =  n
            END IF
  500    CONTINUE

!  check that there is at least one set

         IF ( nsets /= 0 ) THEN

!  print output.

!dir$ ivdep
            DO 510 i = 1, n
               IWK( lsiset + i ) = MIN( IWK( lsiset + i ), nsets + 1 )
!              write( 6, * ) ' variable ', i, ' set ',
!    *                       iwk( lsiset + i )
  510       CONTINUE

! ---------------------------------------------------------------------
!  obtain a list, iwk(lsend), of the variables corresponding to each set
! ---------------------------------------------------------------------

!  clear iwk( lssvse ).

!dir$ ivdep
            DO 520 j = 2, nsets + 2
              IWK( lssvse + j ) = 0
  520       CONTINUE

!  count the number of elements in each set and store in iwk( lssvse ).
!  negate the set numbers in iwk( lsiset ), so that they are flagged
!  as iwk( lsiset ) is gradually overwritten by variable indices

            DO 530 k = 1, n
               j = IWK( lsiset + k )
               IWK( lsiset + k ) = - j
               IWK( lssvse + j + 1 ) = IWK( lssvse + j + 1 ) + 1
  530       CONTINUE

!  compute the starting addresses for each set within iwk(lsiset)

            IWK( lssvse + 1 ) = 1
            DO 540 j = 2, nsets + 2
               IWK( lssvse + j ) = IWK( lssvse + j ) + &
                                   IWK( lssvse + j - 1 )
  540       CONTINUE

!  store in iwk( lsend ) the variable whose set number
!  is the iwk(lssvse + j)-th entry of iwk( lsend )

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

!  reorder the elements into set order.
!  fill in each set from the front. as a new entry is placed
!  in set k increase the pointer IWK(lssvse + k) by one and find
!  the new variable, IWK(lsend + k), that corresponds to the set now
!  pointed to by IWK(lssvse + k)

            DO 660 j = 1, nsets + 1

!  determine the next unplaced entry, IWK(lssvse), in IWK(lsiset)

  610          CONTINUE
               istrt = IWK( lssvse + j )

!  see if all the elements in set j have been assigned

               IF ( istrt == IWK( lssvse + j + 1 ) ) GO TO 660
               IF ( IWK( lsiset + istrt ) > 0 ) GO TO 660

!  extract the variable and set numbers of the starting element

               IVAR = IWK( lsend + j )
               jset = - IWK( lsiset + istrt )

!  move elements in a cycle, ending back at set j

               DO 640 k = istrt, n

!  find the first empty location in set jset in IWK(lsend)

                  inext = IWK( lssvse + jset )

!  extract the variable index of the next element

                  newvar = IWK( lsend + jset )

!  update IWK(lssvse + jset), find the new variable index and store
!  it in iwk(lsend + jset)

                  istarj = inext + 1
                  IWK( lssvse + jset ) = istarj
                  DO 620 ivarp1 = newvar + 1, n
                     IF ( istarj < ivarp1 ) GO TO 630
  620             CONTINUE
                  ivarp1 = n + 1
  630             CONTINUE
                  IWK( lsend + jset ) = ivarp1 - 1

!  if the entry belongs in the j-th set, the cycle is complete

                  IF ( jset == j ) GO TO 650

!  extract the number of the set of the next element

                  newset = - IWK( lsiset + inext )

!  store the variable index of the current element

                  IWK( lsiset + inext ) = IVAR

!  make the next element into the current one

                  IVAR = newvar
                  jset = newset
  640          CONTINUE

!  the cycle is complete

  650          CONTINUE

!  store the variable index of the starting element

               IWK( lsiset + istrt ) = IVAR
               GO TO 610
  660       CONTINUE

!  revise IWK(lssvse) to point to the start of each set

            DO 670 j = nsets + 1, 2, - 1
               IWK( lssvse + j ) = IWK( lssvse + j - 1 )
  670       CONTINUE
            IWK( lssvse + 1 ) = 1
         END IF
!        do 671 i = 1, nsets
!           write( 6, * ) ' set ', i, ' variables ',
!    * ( iwk( lsiset + j ), j = iwk( lssvse + i ),
!    *         iwk( lssvse + i + 1 ) - 1 )
! 671    continue
      ELSE

!  exact gradients are used. no further partitioning of the workspace is needed

         lstype = lsend
         lsswtr = lstype
         lssiwt = lsswtr
         lsiwtr = lssiwt
         lsiset = lsiwtr
         lssvse = lsiset
         lswtra = lwkstr
      END IF

!  set the length of the remaining partitions of the workspace for
!  array bound checking in calls to other subprograms

      lntype = MAX( 1, lsswtr - lstype )
      lnswtr = MAX( 1, lssiwt - lsswtr )
      lnsiwt = MAX( 1, lsiwtr - lssiwt )
      lniwtr = MAX( 1, lsiset - lsiwtr )
      lniset = MAX( 1, lssvse - lsiset )
      lnsvse = MAX( 1, lsend - lssvse )
      lnwtra = MAX( 1, lwkstr - lswtra )

!  record the lengths of the remaining integer and real workspace

      liwk2 = liwk - lsend
      lwk2 = lwk - lwkstr

! -- set the starting addresses for the partitions within fuvals. --

!  a full description of the partitions of fuvals is given in the
!  the introductory comments to subroutine sbmin

      lfxi = 0
      lgxi = lfxi + nel
      lhxi = INTVAR( nel1 ) - 1
      lggfx = lggfx - 1
      ldx = lggfx + n
      lgrjac = ldx + n
      lend = lgrjac + nvargp

!  print all of the starting addresses for the workspace array partitions

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

!  check that the array fuvals has sufficient room for the calculation

      IF ( lend > lfuval ) THEN
         WRITE( iout, 2050 ) lend - lfuval
         status = 5
         RETURN
      END IF

!  set the length of each partition of the real workspace array
!  fuvals for array bound checking in calls to other subprograms

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
      status = 0

!  non-executable statements

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

!  end of subroutine initw

      END SUBROUTINE INITW
