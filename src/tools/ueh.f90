! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE UEH( data, n, X, ne, IRNHI, lirnhi, le, &
                         IPRNHI, HI, lhi, IPRHI, BYROWS )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n, ne, le, lirnhi, lhi 
      LOGICAL :: BYROWS
      INTEGER :: IRNHI ( lirnhi )
      INTEGER :: IPRNHI( le ), IPRHI ( le )
      REAL ( KIND = wp ) :: X ( n ), HI ( lhi )
                         

!  Compute the Hessian matrix of a group partially separable function
!  initially written in Standard Input Format (SIF).

!  The matrix is represented in "finite element format", i.e., 

!           ne
!      H = sum H_i, 
!          i=1

!  where each element H_i involves a small subset of the rows of H.
!  H is stored as a list of the row indices involved in each element
!  and the upper triangle of H_i (stored by rows or columns). 
!  Specifically,

!  ne (integer) number of elements
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

      INTEGER :: i, j, ifstat, igstat
      INTEGER :: ig, liwkh, inform
      REAL ( KIND = wp ) :: zero, one, ftt
      PARAMETER ( zero = 0.0_wp, one = 1.0_wp )
!D    EXTERNAL           RANGE, ELFUN, GROUP, ELGRD, ASMBE

!  there are non-trivial group functions.

      DO 10 i = 1, MAX( data%nel, data%ng )
        data%ICALCF( i ) = i
   10 CONTINUE

!  evaluate the element function values.

      CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), data%nel, &
                   data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                   data%IELVAR( 1 ), data%INTVAR( 1 ), &
                   data%ISTADH( 1 ), data%ISTEP( 1 ), &
                   data%ICALCF( 1 ),  &
                   data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                   data%lntvar, data%lintre, lfuval, data%lvscal, data%lepvlu,  &
                   1, ifstat )

!  evaluate the element function gradients and Hessians.

      CALL ELFUN ( data%FUVALS, X, data%EPVALU( 1 ), data%nel, &
                   data%ITYPEE( 1 ), data%ISTAEV( 1 ), &
                   data%IELVAR( 1 ), data%INTVAR( 1 ), &
                   data%ISTADH( 1 ), data%ISTEP( 1 ), &
                   data%ICALCF( 1 ),  &
                   data%lintre, data%lstaev, data%lelvar, data%lntvar, data%lstadh,  &
                   data%lntvar, data%lintre, lfuval, data%lvscal, data%lepvlu,  &
                   3, ifstat )

!  compute the group argument values ft.

      DO 70 ig = 1, data%ng
         ftt = - data%B( ig )

!  include the contribution from the linear element.

         DO 30 j = data%ISTADA( ig ), data%ISTADA( ig + 1 ) - 1
            ftt = ftt + data%A( j ) * X( data%ICNA( j ) )
   30    CONTINUE

!  include the contributions from the nonlinear elements.

         DO 60 j = data%ISTADG( ig ), data%ISTADG( ig + 1 ) - 1
            ftt = ftt + data%ESCALE( j ) * data%FUVALS( data%IELING( j ) )
   60    CONTINUE
         data%FT( ig ) = ftt

!  Record the derivatives of trivial groups.

         IF ( data%GXEQX( ig ) ) THEN
            data%GVALS( data%ng + ig ) = one
            data%GVALS( 2 * data%ng + ig ) = zero
         END IF
   70 CONTINUE

!  evaluate the group derivative values.

      IF ( .NOT. data%altriv ) CALL GROUP ( data%GVALS( 1 ), data%ng, &
            data%FT( 1 ), data%GPVALU( 1 ), data%ng, &
            data%ITYPEG( 1 ), data%ISTGP( 1 ), &
            data%ICALCF( 1 ), &
            data%lcalcg, data%ng1, data%lcalcg, data%lcalcg, data%lgpvlu, &
            .TRUE., igstat )

!  Compute the gradient value.

      CALL ELGRD( n, data%ng, data%firstg, data%ICNA( 1 ), data%licna, &
                   data%ISTADA( 1 ), data%lstada, data%IELING( 1 ), &
                   data%leling, data%ISTADG( 1 ), data%lstadg, &
                   data%ITYPEE( 1 ), data%lintre, &
                   data%ISTAEV( 1 ), data%lstaev, data%IELVAR( 1 ), &
                   data%lelvar, data%INTVAR( 1 ), data%lntvar, data%IWORK( data%lsvgrp + 1 ), &
                   data%lnvgrp, data%IWORK( data%lstajc + 1 ), data%lnstjc, &
                   data%IWORK( data%lstagv + 1 ), data%lnstgv, data%A( 1 ), data%la, &
                   data%GVALS( data%ng + 1 ), data%lgvals, &
                   data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ), &
                   data%GSCALE( 1 ), data%lgscal, &
                   data%ESCALE( 1 ), data%lescal, data%FUVALS( data%lgrjac + 1 ), &
                   data%lngrjc, data%WRK( 1 ), data%WRK( n + 1 ), data%maxsel, &
                   data%GXEQX( 1 ), data%lgxeqx, &
                   data%INTREP( 1 ), data%lintre, RANGE )
      data%firstg = .FALSE.

!  Define the real work space needed for ASMBE.
!  Ensure that there is sufficient space.

      IF ( data%lwk2 < n + 3 * data%maxsel ) THEN
         IF ( iout > 0 ) WRITE( iout, 2000 )
         STOP
      END IF

!  Define the integer work space needed for ASMBE.
!  Ensure that there is sufficient space.

      liwkh = data%liwk2 - n

!  Assemble the Hessian.

      CALL ASMBE( n, data%ng, data%maxsel,  &
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
                   data%IWORK( liwkh + 1 ), data%liwk2 - liwkh, &
                   data%A( 1 ), data%la, data%FUVALS, data%lnguvl, data%FUVALS, data%lnhuvl, &
                   data%GVALS( data%ng + 1 ), data%GVALS( 2 * data%ng + 1 ), &
                   data%GSCALE( 1 ), data%ESCALE( 1 ), data%lescal, &
                   data%WRK( 1 ), data%lwk2 - data%ng, &
                   data%GXEQX( 1 ), data%lgxeqx, data%INTREP( 1 ), &
                   data%lintre, data%ITYPEE( 1 ), data%lintre, RANGE, ne, &
                   IRNHI, lirnhi, IPRNHI, HI, lhi, IPRHI, &
                   BYROWS, 1, iout, inform )

!  Check that there is room for the elements.

      IF ( inform > 0 ) THEN
         IF ( iout > 0 ) WRITE( iout, 2020 )
         STOP
      END IF

!  Update the counters for the report tool.

      data%nc2oh = data%nc2oh + 1
      RETURN

! Non-executable statements.

 2000 FORMAT( ' ** SUBROUTINE UEH: Increase the size of WK ' )
 2020 FORMAT( ' ** SUBROUTINE UEH: Increase the size of', &
              ' IPNRHI, IPRHI, IRNHI or HI ' )

!  end of UEH.

      END
