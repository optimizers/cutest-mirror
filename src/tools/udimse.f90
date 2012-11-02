! ( Last modified on 10 Sepc 2004 at 16:45:38 )
!  Correction: 10/Sep/2004: undeclared integer variable declared
      SUBROUTINE UDIMSE( data, NE, NZH, NZIRNH )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: NE, NZH, NZIRNH

!  Compute the number of elements and the space required to store
!  the Hessian matrix of the objective function of a problem
!  initially written in Standard Input Format (SIF).

!  The matrix is represented in "finite element format", i.e., 

!           ne
!      H = sum H_i, 
!          i=1

!  where each element H_i involves a small subset of the rows of H.
!  H is stored as a list of the row indices involved in each element
!  and the upper triangle of H_i (stored by rows or columns). 

!  NE (integer) number of elements
!  NZH (integer) number of entries needed to store the real values of
!                   H. Specifically, the sum of the number of entries in
!                   the upper triangle of each H_i.
!  NZIRNH (integer) number of entries needed to store the integer entries
!                   of H. Specifically, the sum of the row dimensions of 
!                   each H_i.

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


!  Integer variables from the PRFCTS common block.


!  Local variables

      INTEGER :: IG, NVARG, IG1

!  Initilaize counts

      NE = 0
      NZH = 0
      NZIRNH = 0

!  Loop over the groups

      DO 10 IG = 1, data%ng
         IG1 = IG + 1

!  Only consider nonlinear groups

         IF ( data%ISTADG( IG ) < data%ISTADG( IG1 ) .OR.  &
              .NOT. data%GXEQX( IG ) ) THEN
            NE = NE + 1
            NVARG = data%IWORK( data%lstagv + IG1 ) - data%IWORK( data%lstagv + IG )
            NZIRNH = NZIRNH + NVARG
            NZH = NZH + ( NVARG * ( NVARG + 1 ) ) / 2
         END IF
   10 CONTINUE
      RETURN

!  end of UDIMSE.

      END
