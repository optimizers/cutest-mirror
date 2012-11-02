! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CSGREH( data, N, M, X, GRLAGF, LV, V,  &
                         NNZJ, LCJAC, CJAC, INDVAR, INDFUN,  &
                         NE, IRNHI, LIRNHI, LE, &
                         IPRNHI, HI, LHI, IPRHI, BYROWS )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: N, M, LV, NNZJ, LCJAC 
      INTEGER :: NE, LE, LIRNHI, LHI 
      LOGICAL :: GRLAGF, BYROWS
      INTEGER :: INDVAR( LCJAC), INDFUN( LCJAC )
      INTEGER :: IRNHI ( LIRNHI )
      INTEGER :: IPRNHI( LE ), IPRHI ( LE )
      REAL ( KIND = wp ) :: X ( N ), V ( LV ),  &
                         HI ( LHI ), CJAC ( LCJAC )

!  Compute the Jacobian matrix for an optimization problem
!  initially written in Standard Input Format (SIF).
!  Also compute the Hessian matrix of the Lagrangian function of
!  the problem.

!  The Jacobian is represented in "co-ordinate" format.
!  The Hessian is represented in "finite element format", i.e., 

!           ne
!      H = sum H_i, 
!          i=1

!  where each element H_i involves a small subset of the rows of H.
!  H is stored as a list of the row indices involved in each element
!  and the upper triangle of H_i (stored by rows or columns). 
!  Specifically,

!  NE (integer) number of elements
!  IRNHI (integer array) a list of the row indices involved which each
!          element. Those for element i directly proceed those for 
!          element i + 1, i = 1, ..., NE-1
!  IPRNHI (integer array) pointers to the position in IRNHI of the first 
!          row index in each element. IPRNHI(NE + 1) points to the first 
!          empty location in IRPNHI
!  HI (real array) a list of the nonzeros in the upper triangle of
!          H_i, stored by rows, or by columns, for each element. Those 
!          for element i directly proceed those for element, i + 1, 
!          i = 1, ..., NE-1
!  IPRHI (integer array) pointers to the position in HI of the first 
!          nonzero in each element. IPRHI(NE + 1) points to the first 
!          empty location in HI
!  BYROWS (logical) must be set .TRUE. if the upper triangle of each H_i 
!          is to be stored by rows, and .FALSE. if it is to be stored
!          by columns.

!  Based on the minimization subroutine data%laNCELOT/SBMIN
!  by Conn, Gould and Toint.

!  Nick Gould, for CGT productions,
!  November 1994.


! ---------------------------------------------------------------------




! ---------------------------------------------------------------------



! ---------------------------------------------------------------------


! ---------------------------------------------------------------------

!  integer variables from the GLOBAL common block.


!  integer variables from the LOCAL common block.


!  integer variables from the PRFCTS common block.


!  Local variables

      INTEGER :: I, J, IEL, K, IG, II, IG1, L, JJ, LL
      INTEGER :: LIWKH, ICON, LGTEMP, INFORM, IENDGV
      INTEGER :: NIN, NVAREL, NELOW, NELUP, ISTRGV
      INTEGER :: IFSTAT, IGSTAT
      LOGICAL :: NONTRV
      REAL ( KIND = wp ) :: FTT, ONE, ZERO, GI, SCALEE, GII
      PARAMETER ( ZERO = 0.0_wp, ONE = 1.0_wp )
      EXTERNAL :: RANGE, ELFUN, GROUP 
!D    EXTERNAL           DSETVL, DSETVI, DELGRD, DASMBE

!  there are non-trivial group functions.

      DO 10 I = 1, MAX( data%nelnum, data%ng )
        data%ICALCF( I ) = I
   10 CONTINUE

!  evaluate the element function values.

      CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), data%nelnum, &
                   data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                   data%IELVAR( 1 ), data%INTVAR( 1 ), &
                   data%ISTADH( 1 ), data%ISTEP( 1 ), &
                   data%ICALCF( 1 ),  &
                   data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                   data%lntvar, data%lintre, LFUVAL, data%lvscal, data%lepvlu,  &
                   1, IFSTAT )

!  evaluate the element function gradients and Hessians.

      CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), data%nelnum, &
                   data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                   data%IELVAR( 1 ), data%INTVAR( 1 ), &
                   data%ISTADH( 1 ), data%ISTEP( 1 ), &
                   data%ICALCF( 1 ),  &
                   data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                   data%lntvar, data%lintre, LFUVAL, data%lvscal, data%lepvlu,  &
                   3, IFSTAT )

!  compute the group argument values ft.

      DO 70 IG = 1, data%ng
         FTT = - data%B( IG )

!  include the contribution from the linear element.

         DO 30 J = data%ISTADA( IG ), data%ISTADA( IG + 1 ) - 1
            FTT = FTT + data%A( J ) * X( data%ICNA( J ) )
   30    CONTINUE

!  include the contributions from the nonlinear elements.

         DO 60 J = data%ISTADG( IG ), data%ISTADG( IG + 1 ) - 1
            FTT = FTT + data%ESCALE( J ) * data%FUVALS( data%IELING( J ) )
   60    CONTINUE
         data%FT( IG ) = FTT

!  Record the derivatives of trivial groups.

         IF ( data%GXEQX( IG ) ) THEN
            data%GVALS( data%ng + IG ) = ONE
            data%GVALS( 2 * data%ng + IG ) = ZERO
         END IF
   70 CONTINUE

!  evaluate the group derivative values.

      IF ( .NOT. data%altriv ) CALL GROUP ( data%GVALS( 1 ), data%ng, &
            data%FT( 1 ), data%GPVALU( 1 ), data%ng, &
            data%ITYPEG( 1 ), data%ISTGP( 1 ), &
            data%ICALCF( 1 ), &
            data%lcalcg, data%ng1, data%lcalcg, data%lcalcg, data%lgpvlu, &
            .TRUE., IGSTAT )

!  Define the real work space needed for ELGRD.
!  Ensure that there is sufficient space.

      IF ( data%lwk2 < data%ng ) THEN
         IF ( IOUT > 0 ) WRITE( IOUT, 2000 )
         STOP
      END IF
      IF ( data%numcon > 0 ) THEN

!  Change the group weightings to include the contributions from
!  the Lagrange multipliers.

         DO 80 IG = 1, data%ng
            I = data%KNDOFC( IG )
            IF ( I == 0 ) THEN
               data%WRK( IG ) = data%GSCALE( IG )
            ELSE
               data%WRK( IG ) = data%GSCALE( IG ) * V( I )
            END IF
   80    CONTINUE

!  Compute the gradient values. Initialize the gradient of the
!  objective function as zero.

         NNZJ = 0
         LGTEMP = WRK + N + data%maxsel
         DO 120 J = 1, N
            data%WRK( LGTEMP + J ) = ZERO
  120    CONTINUE

!  Consider the IG-th group.

         DO 290 IG = 1, data%ng
            IG1 = IG + 1
            ICON = data%KNDOFC( IG )
            ISTRGV = data%IWORK( data%lstagv + IG )
            IENDGV = data%IWORK( data%lstagv + IG1 ) - 1
            NELOW = data%ISTADG( IG )
            NELUP = data%ISTADG( IG1 ) - 1
            NONTRV = .NOT. data%GXEQX( IG )

!  Compute the first derivative of the group.

            GI = data%GSCALE( IG )
            GII = data%WRK( IG )
            IF ( NONTRV ) THEN
               GI = GI  * data%GVALS( data%ng + IG )
               GII = GII * data%GVALS( data%ng + IG )
            END IF
      CALL DSETVI( IENDGV - ISTRGV + 1, data%WRK( 1 ), &
                         data%IWORK( data%lsvgrp + ISTRGV ), ZERO )

!  This is the first gradient evaluation or the group has nonlinear
!  elements.

            IF ( data%firstg .OR. NELOW <= NELUP ) THEN

!  Loop over the group's nonlinear elements.

               DO 150 II = NELOW, NELUP
                  IEL = data%IELING( II )
                  K = data%INTVAR( IEL )
                  L = data%ISTAEV( IEL )
                  NVAREL = data%ISTAEV( IEL + 1 ) - L
                  SCALEE = data%ESCALE( II )
                  IF ( data%INTREP( IEL ) ) THEN

!  The IEL-th element has an internal representation.

                     NIN = data%INTVAR( IEL + 1 ) - K
                     CALL RANGE ( IEL, .TRUE., data%FUVALS( K ), &
                                  data%WRK( N + 1 ), NVAREL, NIN, &
                                  data%ITYPEE( IEL ), &
                                  NIN, NVAREL )
!DIR$ IVDEP
                     DO 130 I = 1, NVAREL
                        J = data%IELVAR( L )
                        data%WRK( J ) = data%WRK( J ) + &
                                           SCALEE * data%WRK( N + I )
                        L = L + 1
  130                CONTINUE
                  ELSE

!  The IEL-th element has no internal representation.

!DIR$ IVDEP
                     DO 140 I = 1, NVAREL
                        J = data%IELVAR( L )
                        data%WRK( J ) = data%WRK( J ) + &
                                           SCALEE * data%FUVALS( K )
                        K = K + 1
                        L = L + 1
  140                CONTINUE
                  END IF
  150          CONTINUE

!  Include the contribution from the linear element.

!DIR$ IVDEP
               DO 160 K = data%ISTADA( IG ), &
                                  data%ISTADA( IG1 ) - 1
                  J = data%ICNA( K )
                  data%WRK( J ) = data%WRK( J ) + data%A( K )
  160          CONTINUE

!  Allocate a gradient.

!DIR$ IVDEP
               DO 190 I = ISTRGV, IENDGV
                  LL = data%IWORK( data%lsvgrp + I )

!  The group belongs to the objective function.

                  IF ( ICON == 0 ) THEN
                     data%WRK( LGTEMP + LL ) = data%WRK( LGTEMP + LL ) + &
                                         GI * data%WRK( LL )

!  The group defines a constraint.

                  ELSE
                     NNZJ = NNZJ + 1
                     IF ( NNZJ <= LCJAC ) THEN
                        CJAC ( NNZJ ) = GI * data%WRK( LL )
                        INDFUN( NNZJ ) = ICON
                        INDVAR( NNZJ ) = LL
                     END IF
                     IF ( GRLAGF ) &
                        data%WRK( LGTEMP + LL ) = data%WRK( LGTEMP + LL ) + &
                                            GII * data%WRK( LL )
                  END IF

!  If the group is non-trivial, also store the nonzero entries of the
!  gradient of the function in GRJAC.

                  IF ( NONTRV ) THEN
                     JJ = data%IWORK( data%lstajc + LL )
                     data%FUVALS( data%lgrjac + JJ ) = data%WRK( LL )

!  Increment the address for the next nonzero in the column of
!  the jacobian for variable LL.

                     data%IWORK( data%lstajc + LL ) = JJ + 1
                  END IF
  190          CONTINUE

!  This is not the first gradient evaluation and there is only a linear
!  element.

            ELSE
!                            linear element improved. 43 lines replace 40

!  Include the contribution from the linear element.

!DIR$ IVDEP
               DO 210 K = data%ISTADA( IG ),data%ISTADA( IG1 ) - 1
                  J = data%ICNA( K )
                  data%WRK( J ) = data%WRK( J ) + data%A( K )
  210          CONTINUE

!  Allocate a gradient.

!DIR$ IVDEP
               DO 220 I = ISTRGV, IENDGV
                  LL = data%IWORK( data%lsvgrp + I )

!  The group belongs to the objective function.

                  IF ( ICON == 0 ) THEN
                     data%WRK( LGTEMP + LL ) = data%WRK( LGTEMP + LL ) + &
                                         GI * data%WRK( LL )

!  The group defines a constraint.

                  ELSE
                     NNZJ = NNZJ + 1
                     IF ( NNZJ <= LCJAC ) THEN
                        CJAC ( NNZJ ) = GI * data%WRK( LL )
                        INDFUN( NNZJ ) = ICON
                        INDVAR( NNZJ ) = LL
                     END IF
                     IF ( GRLAGF ) &
                        data%WRK( LGTEMP + LL ) = data%WRK( LGTEMP + LL ) + &
                                            GII * data%WRK( LL )
                  END IF

!  Increment the address for the next nonzero in the column of
!  the jacobian for variable LL.

                  IF ( NONTRV ) THEN
                     JJ = data%IWORK( data%lstajc + LL )
                     data%IWORK( data%lstajc + LL ) = JJ + 1
                  END IF
  220          CONTINUE
            END IF
  290    CONTINUE

!  Reset the starting addresses for the lists of groups using
!  each variable to their values on entry.

         DO 300 I = N, 2, - 1
            data%IWORK( data%lstajc + I ) = data%IWORK( data%lstajc + I - 1 )
  300    CONTINUE
         data%IWORK( data%lstajc + 1 ) = 1

!  Transfer the gradient of the objective function to the sparse
!  storage scheme.

         DO 310 I = 1, N
!           IF ( data%WRK( LGTEMP + I ) /= ZERO ) THEN
               NNZJ = NNZJ + 1
               IF ( NNZJ <= LCJAC ) THEN
                  CJAC ( NNZJ ) = data%WRK( LGTEMP + I )
                  INDFUN( NNZJ ) = 0
                  INDVAR( NNZJ ) = I
               END IF
!           END IF
  310    CONTINUE
      ELSE

!  Compute the gradient value.

      CALL DELGRD( N, data%ng, data%firstg, data%ICNA( 1 ), data%licna, &
                      data%ISTADA( 1 ), data%lstada, data%IELING( 1 ), &
                      data%leling, data%ISTADG( 1 ), data%lstadg, &
                      data%ITYPEE( 1 ), data%lintre, &
                      data%ISTAEV( 1 ), data%lstaev, data%IELVAR( 1 ), &
                      data%lelvar, data%INTVAR( 1 ), data%lntvar, &
                      data%IWORK( data%lsvgrp + 1 ), &
                      data%lnvgrp, data%IWORK( data%lstajc + 1 ), data%lnstjc, &
                      data%IWORK( data%lstagv + 1 ), data%lnstgv, data%A( 1 ), data%la, &
                      data%GVALS( data%ng + 1 ), data%lgvals, &
                      data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ), &
                      data%GSCALE( 1 ), data%lgscal, &
                      data%ESCALE( 1 ), data%lescal, data%FUVALS( data%lgrjac + 1 ), &
                      data%lngrjc, data%WRK( 1 ), data%WRK( N + 1 ), data%maxsel, &
                      data%GXEQX( 1 ), data%lgxeqx, &
                      data%INTREP( 1 ), data%lintre, RANGE )

!  Transfer the gradient of the objective function to the sparse
!  storage scheme.

         NNZJ = 0
         DO 400 I = 1, N
!           IF ( data%FUVALS( data%lggfx + I ) /= ZERO ) THEN
               NNZJ = NNZJ + 1
               IF ( NNZJ <= LCJAC ) THEN
                  CJAC ( NNZJ ) = data%FUVALS( data%lggfx + I )
                  INDFUN( NNZJ ) = 0
                  INDVAR( NNZJ ) = I
               END IF
!           END IF
  400    CONTINUE
      END IF
      data%firstg = .FALSE.

!  Verify that the Jacobian can fit in the alloted space

      IF ( NNZJ > LCJAC ) THEN
         IF ( IOUT > 0 ) WRITE( IOUT, 2030 ) NNZJ - LCJAC 
         STOP
      END IF

!  Define the real work space needed for ASMBE.
!  Ensure that there is sufficient space.

      IF ( data%numcon > 0 ) THEN
         IF ( data%lwk2 < N + 3 * data%maxsel + data%ng ) THEN
            IF ( IOUT > 0 ) WRITE( IOUT, 2000 )
            STOP
         END IF
      ELSE
         IF ( data%lwk2 < N + 3 * data%maxsel ) THEN
            IF ( IOUT > 0 ) WRITE( IOUT, 2000 )
            STOP
         END IF
      END IF

!  Define the integer work space needed for ASMBE.
!  Ensure that there is sufficient space.

      LIWKH = data%liwk2 - N

!  Assemble the Hessian.

      IF ( data%numcon > 0 ) THEN
      CALL DASMBE( N, data%ng, data%maxsel,  &
                      data%ISTADH( 1 ), data%lstadh, &
                      data%ICNA( 1 ), data%licna, &
                      data%ISTADA( 1 ), data%lstada, &
                      data%INTVAR( 1 ), data%lntvar, &
                      data%IELVAR( 1 ), data%lelvar, &
                      data%IELING( 1 ), data%leling, &
                      data%ISTADG( 1 ), data%lstadg, &
                      data%ISTAEV( 1 ), data%lstaev, &
                      data%IWORK( data%lstagv + 1 ), data%lnstgv, &
                      data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, &
                      data%IWORK( LIWKH + 1 ), data%liwk2 - LIWKH, &
                      data%A( 1 ), data%la, data%FUVALS, data%lnguvl, data%FUVALS, data%lnhuvl, &
                      data%GVALS( data%ng + 1 ), data%GVALS( 2 * data%ng + 1 ), &
                      data%WRK( 1 ), data%ESCALE( 1 ), data%lescal, &
                      data%WRK( data%ng + 1 ), data%lwk2 - data%ng, &
                      data%GXEQX( 1 ), data%lgxeqx, data%INTREP( 1 ), &
                      data%lintre, data%ITYPEE( 1 ), data%lintre, RANGE, NE,  &
                      IRNHI, LIRNHI, IPRNHI, HI, LHI, IPRHI, &
                      BYROWS, 1, IOUT, INFORM )
      ELSE
      CALL DASMBE( N, data%ng, data%maxsel,  &
                      data%ISTADH( 1 ), data%lstadh, &
                      data%ICNA( 1 ), data%licna, &
                      data%ISTADA( 1 ), data%lstada, &
                      data%INTVAR( 1 ), data%lntvar, &
                      data%IELVAR( 1 ), data%lelvar, &
                      data%IELING( 1 ), data%leling, &
                      data%ISTADG( 1 ), data%lstadg, &
                      data%ISTAEV( 1 ), data%lstaev, &
                      data%IWORK( data%lstagv + 1 ), data%lnstgv, &
                      data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, &
                      data%IWORK( LIWKH + 1 ), data%liwk2 - LIWKH, &
                      data%A( 1 ), data%la, data%FUVALS, data%lnguvl, data%FUVALS, data%lnhuvl, &
                      data%GVALS( data%ng + 1 ), data%GVALS( 2 * data%ng + 1 ), &
                      data%GSCALE( 1 ), data%ESCALE( 1 ), data%lescal, &
                      data%WRK( 1 ), data%lwk2 - data%ng, &
                      data%GXEQX( 1 ), data%lgxeqx, data%INTREP( 1 ), &
                      data%lintre, data%ITYPEE( 1 ), data%lintre, RANGE, NE,  &
                      IRNHI, LIRNHI, IPRNHI, HI, LHI, IPRHI, &
                      BYROWS, 1, IOUT, INFORM )
      END IF

!  Check that there is room for the elements

      IF ( INFORM > 0 ) THEN
         IF ( IOUT > 0 ) WRITE( IOUT, 2020 )
         STOP
      END IF

!  Update the counters for the report tool.

      data%nc2cg = data%nc2og + 1
      data%nc2oh = data%nc2oh + 1
      data%nc2cg = data%nc2cg + data%pnc
      data%nc2ch = data%nc2ch + data%pnc
      RETURN

! Non-executable statements.

 2000 FORMAT( ' ** SUBROUTINE CSGREH: Increase the size of WK ' )
 2020 FORMAT( ' ** SUBROUTINE CSGREH: Increase the size of', &
              ' IPNRHI, IPRHI, IRNHI or HI ' )
 2030 FORMAT( /, ' ** SUBROUTINE CSGREH: array length LCJAC too small.', &
              /, ' -- Minimization abandoned.', &
              /, ' -- Increase the parameter LCJAC by at least ', I8, &
                 ' and restart.' )

!  end of CSGREH.

      END
