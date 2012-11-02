! ( Last modified on 10 Sepc 2004 at 16:45:38 )
!  Correction: 10/Sep/2004: undeclared integer variables declared
      SUBROUTINE CIDH ( data, N, X, IPROB, LH1, H )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: N, IPROB, LH1
      REAL ( KIND = wp ) :: X ( N ), H( LH1, N )

!  Compute the Hessian matrix of a specified problem function 
! (IPROB = 0 is the objective function, while IPROB > 0 is the 
!   IPROB-th constraint) of a problem initially written in 
!  Standard Input Format (SIF).

!  H is a two-dimensional array which gives the value of the
!    Hessian matrix of the specified problem function evaluated at
!    X. The i,j-th component of the array will contain
!    the derivative with respect to variables X(i) and X(j).

!  Based on the minimization subroutine data%laNCELOT/SBMIN
!  by Conn, Gould and Toint.

!  Nick Gould, for CGT productions,
!  August 1998.


! ---------------------------------------------------------------------




! ---------------------------------------------------------------------



! ---------------------------------------------------------------------


! ---------------------------------------------------------------------

!  integer variables from the GLOBAL common block.


!  integer variables from the LOCAL common block.


!  integer variables from the PRFCTS common block.


!  Local variables

      INTEGER :: I, J, NCALCF, NCALCG, IFSTAT, IGSTAT
      INTEGER :: IG, LIH, LIWKH, LNXTRW, LINXTR, INFORM
      INTEGER :: K, LH, LIRNH, LJCNH, NNZH, LWKH
      REAL ( KIND = wp ) :: ZERO, ONE, FTT
      PARAMETER ( ZERO = 0.0_wp, ONE = 1.0_wp )
      EXTERNAL :: RANGE 

!  Check input parameters.

      IF ( IPROB < 0 ) THEN
         WRITE( IOUT, 2020 )
         STOP
      END IF

!  Find group index IG of constraint IPROB.

      IF ( IPROB > 0 ) THEN
         IG = 0
         DO 10 I = 1, data%ng
            IF ( data%KNDOFC( I ) == IPROB ) THEN
              IG = I
              GO TO 20
            END IF
   10    CONTINUE
   20    CONTINUE
         IF ( IG == 0 ) THEN
            WRITE( IOUT, 2020 )
            STOP
         END IF
      END IF

!  Find which elements are involved in the required problem function.
!  Initialize the list to zero.

      DO 30 I = 1, data%nelnum
        data%ICALCF( I ) = 0
   30 CONTINUE

!  If group IG is involved, record its elements.

      DO 50 IG = 1, data%ng
         IF ( data%KNDOFC( IG ) == IPROB ) THEN
            DO 40 J = data%ISTADG( IG ), data%ISTADG( IG + 1 ) - 1
               data%ICALCF( data%IELING( J ) ) = 1
   40       CONTINUE
         END IF
   50 CONTINUE

!  Only compute the elements which are involved in the required function.

      NCALCF = 0
      DO 80 I = 1, data%nelnum
         IF ( data%ICALCF( I ) == 1 ) THEN
            NCALCF = NCALCF + 1
            data%ICALCF( NCALCF ) = I
         ELSE

!  If this is the first ever evaluation, initialize data%FUVALS.

            IF ( data%firstg ) THEN
               data%FUVALS( I ) = ZERO
               DO 60 J = data%INTVAR( I ), data%INTVAR( I + 1 ) - 1
                  data%FUVALS( J ) = ZERO
   60          CONTINUE
               DO 70 J = data%ISTADH( I ), data%ISTADH( I + 1 ) - 1
                  data%FUVALS( J ) = ZERO
   70          CONTINUE
            END IF
         END IF
   80 CONTINUE

!  evaluate the element function values.

      CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), NCALCF, &
                   data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                   data%IELVAR( 1 ), data%INTVAR( 1 ), &
                   data%ISTADH( 1 ), data%ISTEP( 1 ), &
                   data%ICALCF( 1 ),  &
                   data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                   data%lntvar, data%lintre, LFUVAL, data%lvscal, data%lepvlu,  &
                   1, IFSTAT )

!  evaluate the element function derivatives.

      CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), NCALCF, &
                   data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                   data%IELVAR( 1 ), data%INTVAR( 1 ), &
                   data%ISTADH( 1 ), data%ISTEP( 1 ), &
                   data%ICALCF( 1 ),  &
                   data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                   data%lntvar, data%lintre, LFUVAL, data%lvscal, data%lepvlu,  &
                   3, IFSTAT )

!  Compute the list of groups involved in the required problem function.

      NCALCG = 0
      DO 130 IG = 1, data%ng
         IF ( data%KNDOFC( IG ) == IPROB ) THEN
            NCALCG = NCALCG + 1
            data%ICALCF( NCALCG ) = IG

!  compute the group argument values ft.

            FTT = - data%B( IG )

!  include the contribution from the linear element.

            DO 110 J = data%ISTADA( IG ), data%ISTADA( IG + 1 ) - 1
               FTT = FTT + data%A( J ) * X( data%ICNA( J ) )
  110       CONTINUE

!  include the contributions from the nonlinear elements.

            DO 120 J = data%ISTADG( IG ), data%ISTADG( IG + 1 ) - 1
               FTT = FTT + data%ESCALE( J ) *  &
                      data%FUVALS( data%IELING( J ) )
  120          CONTINUE
            data%FT( IG ) = FTT

!  Record the derivatives of trivial groups.

            IF ( data%GXEQX( IG ) ) THEN
               data%GVALS( data%ng + IG ) = ONE
               data%GVALS( 2 * data%ng + IG ) = ZERO
            END IF
         ELSE

!  If this is the first ever evaluation, initialize GVALS.

            IF ( data%firstg ) THEN
               data%GVALS( data%ng + IG ) = ONE
               data%GVALS( 2 * data%ng + IG ) = ZERO
            END IF
         END IF
  130 CONTINUE

!  Evaluate the group derivative values.

      IF ( .NOT. data%altriv ) CALL GROUP ( data%GVALS( 1 ), data%ng, &
            data%FT( 1 ), data%GPVALU( 1 ), NCALCG, &
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

         DO 140 IG = 1, data%ng
            I = data%KNDOFC( IG )
            IF ( I == IPROB ) THEN
               data%WRK( IG ) = data%GSCALE( IG )
            ELSE
               data%WRK( IG ) = ZERO
            END IF
  140    CONTINUE

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
                      data%WRK( 1 ), data%ng, &
                      data%ESCALE( 1 ), data%lescal, data%FUVALS( data%lgrjac + 1 ), &
                      data%lngrjc, data%WRK( 1 ), data%WRK( N + 1 ), data%maxsel, &
                      data%GXEQX( 1 ), data%lgxeqx, &
                      data%INTREP( 1 ), data%lintre, RANGE )
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
      END IF
      data%firstg = .FALSE.

!  Define the real work space needed for ASMBL.
!  Ensure that there is sufficient space.

      IF ( data%numcon > 0 ) THEN
         LWKH = data%lwk2 - N - 3 * data%maxsel - data%ng
      ELSE
         LWKH = data%lwk2 - N - 3 * data%maxsel
      END IF
      IF ( LWKH <= 0 ) THEN
         IF ( IOUT > 0 ) WRITE( IOUT, 2000 )
         STOP
      END IF

!  Define the integer work space needed for ASMBL.
!  Ensure that there is sufficient space.

      LIWKH = data%liwk2 - N
      LH = MIN( LWKH, ( LIWKH - 3 * N ) / 4 )
      LINXTR = LH + N
      IF ( LH <= 0 ) THEN
         IF ( IOUT > 0 ) WRITE( IOUT, 2010 )
         STOP
      END IF

!  Set starting addresses for partitions of the integer workspace.

      LIH = LH
      LIRNH = 0
      LJCNH = LIRNH + LIH
      LNXTRW = LJCNH + LIH
      DO 90 I = 1, N
         data%IVAR( I ) = I
   90 CONTINUE

!  Assemble the Hessian.

      IF ( data%numcon > 0 ) THEN
      CALL DASMBL( N, data%ng, data%maxsel, N, LH, LIH, NNZH, &
                   N, data%IVAR( 1), data%ISTADH( 1 ), data%lstadh, &
                   data%ICNA( 1 ), data%licna, &
                   data%ISTADA( 1 ), data%lstada, data%INTVAR( 1 ), data%lntvar, &
                   data%IELVAR( 1 ), data%lelvar, data%IELING( 1 ), data%leling, &
                   data%ISTADG( 1 ), data%lstadg, data%ISTAEV( 1 ), data%lstaev, &
                   data%IWORK( data%lstagv + 1 ), data%lnstgv, data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, &
                   data%IWORK( data%lsend + LIRNH + 1 ), data%IWORK( data%lsend + LJCNH + 1 ), &
                   data%IWORK( data%lsend + LNXTRW + 1 ), LINXTR, &
                   data%IWORK( data%lsend + LIWKH + 1 ), N, &
                   data%A( 1 ), data%la, data%FUVALS, data%lnguvl, data%FUVALS, data%lnhuvl, &
                   data%GVALS( data%ng + 1 ), data%GVALS( 2 * data%ng + 1 ), &
                   data%WRK( 1 ), data%ESCALE( 1 ), data%lescal, &
                   data%WRK( data%ng + 1 ), data%WRK( LWKH + data%ng + 1 ), &
                   data%lwk2 - LWKH, data%GXEQX( 1 ), &
                   data%lgxeqx, data%INTREP( 1 ), data%lintre, &
                   data%ITYPEE( 1 ), data%lintre, &
                   RANGE, 1, IOUT, .FALSE., I, INFORM, .FALSE., &
                   .FALSE. )
      ELSE
      CALL DASMBL( N, data%ng, data%maxsel, N, LH, LIH, NNZH, &
                   N, data%IVAR( 1), data%ISTADH( 1 ), data%lstadh, &
                   data%ICNA( 1 ), data%licna, &
                   data%ISTADA( 1 ), data%lstada, data%INTVAR( 1 ), data%lntvar, &
                   data%IELVAR( 1 ), data%lelvar, data%IELING( 1 ), data%leling, &
                   data%ISTADG( 1 ), data%lstadg, data%ISTAEV( 1 ), data%lstaev, &
                   data%IWORK( data%lstagv + 1 ), data%lnstgv, data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, &
                   data%IWORK( data%lsend + LIRNH + 1 ), data%IWORK( data%lsend + LJCNH + 1 ), &
                   data%IWORK( data%lsend + LNXTRW + 1 ), LINXTR, &
                   data%IWORK( data%lsend + LIWKH + 1 ), N, &
                   data%A( 1 ), data%la, data%FUVALS, data%lnguvl, data%FUVALS, data%lnhuvl, &
                   data%GVALS( data%ng + 1 ), data%GVALS( 2 * data%ng + 1 ), &
                   data%GSCALE( 1 ), data%ESCALE( 1 ), data%lescal, &
                   data%WRK( 1 ), data%WRK( LWKH + 1 ), &
                   data%lwk2 - LWKH, data%GXEQX( 1 ), &
                   data%lgxeqx, data%INTREP( 1 ), data%lintre, &
                   data%ITYPEE( 1 ), data%lintre, &
                   RANGE, 1, IOUT, .FALSE., I, INFORM, .FALSE., &
                   .FALSE. )
      END IF

!  Check that there is sufficient integer workspace.

      IF ( INFORM > 0 ) THEN
         IF ( IOUT > 0 ) WRITE( IOUT, 2010 )
         STOP
      END IF

!  Initialize the dense matrix.

      DO 220 J = 1, N
         DO 210 I = 1, N
            H( I, J ) = ZERO
  210    CONTINUE
  220 CONTINUE

!  Transfer the matrix from co-ordinate to dense storage and
!  symmetrize the martix.

      IF ( data%numcon > 0 ) THEN
         DO 230 K = 1, NNZH
            I = data%IWORK( data%lsend + LIRNH + K )
            J = data%IWORK( data%lsend + LJCNH + K )
            H( I, J ) = data%WRK( data%ng + K )
            H( J, I ) = data%WRK( data%ng + K )
  230    CONTINUE
      ELSE
         DO 240 K = 1, NNZH
            I = data%IWORK( data%lsend + LIRNH + K )
            J = data%IWORK( data%lsend + LJCNH + K )
            H( I, J ) = data%WRK( K )
            H( J, I ) = data%WRK( K )
  240    CONTINUE
      END IF

!  Update the counters for the report tool.

      IF( IPROB == 0 ) THEN
         data%nc2oh = data%nc2oh + 1
      ELSE 
         data%nc2ch = data%nc2ch + 1
      ENDIF 
      RETURN

! Non-executable statements.

 2000 FORMAT( ' ** SUBROUTINE CIDH: Increase the size of WK ' )
 2010 FORMAT( ' ** SUBROUTINE CIDH: Increase the size of IWK ' )
 2020 FORMAT( ' ** SUBROUTINE CIDH: invalid problem index IPROB ' )

!  end of CIDH.

      END
