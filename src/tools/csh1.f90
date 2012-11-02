! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CSH1 ( data, N, M, X, LV, V, NNZH, &
                         LH, H, IRNH, ICNH )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: N, M, LV, NNZH, LH
      INTEGER :: IRNH ( LH ), ICNH ( LH )
      REAL ( KIND = wp ) :: X ( N ), V ( LV ), H( LH )

!  Compute the Hessian matrix of the constraint part of the Lagrangian 
!  function of a problem initially written in Standard Input Format (SIF).

!  The upper triangle of the Hessian is stored in coordinate form,
!  i.e., the entry H(i) has row index IRNH(i) and column index ICNH(i)
!  for i = 1, ...., NNZH.

!  Based on the minimization subroutine data%laNCELOT/SBMIN
!  by Conn, Gould and Toint.

!  Version which ignores objective function for FilterSQP

!  Nick Gould, for CGT productions,
!  November 1991.


! ---------------------------------------------------------------------




! ---------------------------------------------------------------------


! ---------------------------------------------------------------------


! ---------------------------------------------------------------------

!  integer variables from the GLOBAL common block.


!  integer variables from the LOCAL common block.


!  integer variables from the PRFCTS common block.


!  Local variables

      INTEGER :: I, J, IFSTAT, IGSTAT
      INTEGER :: IG, LIH, LIWKH, LNXTRW, LINXTR, INFORM
      REAL ( KIND = wp ) :: ZERO, ONE, FTT
      PARAMETER ( ZERO = 0.0_wp, ONE = 1.0_wp )
      EXTERNAL :: RANGE 

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

!  evaluate the element function values.

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
!                            for unconstrained problems.
      IF ( data%numcon > 0 ) THEN

!  Change the group weightings to include the contributions from
!  the Lagrange multipliers.

         DO 80 IG = 1, data%ng
            I = data%KNDOFC( IG )
            IF ( I == 0 ) THEN
               data%WRK( IG ) = ZERO
            ELSE
               data%WRK( IG ) = data%GSCALE( IG ) * V( I )
            END IF
   80    CONTINUE

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

!                            for unconstrained problems.
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

!  Define the integer work space needed for ASMBL.
!  Ensure that there is sufficient space.

      LIWKH = data%liwk2 - N
      LINXTR = LIWKH / 2
      IF ( LINXTR < N ) THEN
         IF ( IOUT > 0 ) WRITE( IOUT, 2010 )
         STOP
      END IF

!  Set starting addresses for partitions of the integer workspace.

      LIH = LH
      LNXTRW = 0
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
                   IRNH, ICNH, data%IWORK( data%lsend + LNXTRW + 1 ), LINXTR, &
                   data%IWORK( data%lsend + LIWKH + 1 ), N, &
                   data%A( 1 ), data%la, data%FUVALS, data%lnguvl, data%FUVALS, data%lnhuvl, &
                   data%GVALS( data%ng + 1 ), data%GVALS( 2 * data%ng + 1 ), &
                   data%WRK( 1 ), data%ESCALE( 1 ), data%lescal, &
                   H, data%WRK( data%ng + 1 ), data%lwk2 - data%ng, &
                   data%GXEQX( 1 ), data%lgxeqx, data%INTREP( 1 ), &
                   data%lintre, data%ITYPEE( 1 ), data%lintre, &
                   RANGE, 1, IOUT, .FALSE., I, INFORM,  &
                   .FALSE., .TRUE. )
      ELSE
      CALL DASMBL( N, data%ng, data%maxsel, N, LH, LIH, NNZH, &
                   N, data%IVAR( 1), data%ISTADH( 1 ), data%lstadh, &
                   data%ICNA( 1 ), data%licna, &
                   data%ISTADA( 1 ), data%lstada, data%INTVAR( 1 ), data%lntvar, &
                   data%IELVAR( 1 ), data%lelvar, data%IELING( 1 ), data%leling, &
                   data%ISTADG( 1 ), data%lstadg, data%ISTAEV( 1 ), data%lstaev, &
                   data%IWORK( data%lstagv + 1 ), data%lnstgv, data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, &
                   IRNH, ICNH, data%IWORK( data%lsend + LNXTRW + 1 ), LINXTR, &
                   data%IWORK( data%lsend + LIWKH + 1 ), N,  &
                   data%A( 1 ), data%la, data%FUVALS, data%lnguvl, data%FUVALS, data%lnhuvl, &
                   data%GVALS( data%ng + 1 ), data%GVALS( 2 * data%ng + 1 ), &
                   data%GSCALE( 1 ), data%ESCALE( 1 ), data%lescal, &
                   H, data%WRK( 1 ), data%lwk2 - data%ng, &
                   data%GXEQX( 1 ), data%lgxeqx, data%INTREP( 1 ), &
                   data%lintre, data%ITYPEE( 1 ), data%lintre, &
                   RANGE, 1, IOUT, .FALSE., I, INFORM, &
                   .FALSE., .TRUE. )
      END IF

!  Check that there is sufficient integer workspace.

      IF ( INFORM > 0 ) THEN
         IF ( IOUT > 0 ) WRITE( IOUT, 2010 )
         STOP
      END IF

!  Update the counters for the report tool.

      data%nc2oh = data%nc2oh + 1
      data%nc2ch = data%nc2ch + data%pnc
      RETURN

! Non-executable statements.

 2000 FORMAT( ' ** SUBROUTINE CSH1: Increase the size of WK ' )
 2010 FORMAT( ' ** SUBROUTINE CSH1: Increase the sizes of IWK and LH' )

!  end of CSH1.

      END
