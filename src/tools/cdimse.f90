! ( Last modified on 10 Sepc 2004 at 16:35:38 )
!  Correction: 10/Sep/2004: undeclared integer variables declared
      SUBROUTINE CDIMSE( data, status, ne, nzh, nzirnh )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: ne, nzh, nzirnh
      INTEGER, INTENT( OUT ) :: status

!  Compute the number of elements and the space required to store
!  the Hessian matrix of the Lagrangian function of a problem
!  initially written in Standard Input Format (SIF).

!  The matrix is represented in "finite element format", i.e., 

!           ne
!      H = sum H_i, 
!          i=1

!  where each element H_i involves a small subset of the rows of H.
!  H is stored as a list of the row indices involved in each element
!  and the upper triangle of H_i (stored by rows or columns). 

!  ne (integer) number of elements
!  nzh (integer) number of entries needed to store the real values of
!                   H. Specifically, the sum of the number of entries in
!                   the upper triangle of each H_i.
!  nzirnh (integer) number of entries needed to store the integer entries
!                   of H. Specifically, the sum of the row dimensions of 
!                   each H_i.

!  Based on the minimization subroutine data%laNCELOT/SBMIN
!  by Conn, Gould and Toint.

!  Nick Gould, for CGT productions,
!  November 1994.

!  Local variables

      INTEGER :: ig, nvarg, ig1

!  Initilaize counts

      ne = 0
      nzh = 0
      nzirnh = 0

!  Loop over the groups

      DO ig = 1, data%ng
        ig1 = ig + 1

!  Only consider nonlinear groups

        IF ( data%ISTADG( ig ) < data%ISTADG( ig1 ) .OR.                       &
              .NOT. data%GXEQX( ig ) ) THEN
          ne = ne + 1
          nvarg                                                                &
            = data%IWORK( data%lstagv + ig1 ) - data%IWORK( data%lstagv + ig )
          nzirnh = nzirnh + nvarg
          nzh = nzh + ( nvarg * ( nvarg + 1 ) ) / 2
        END IF
      END DO
      status = 0
      RETURN

!  end of CDIMSE

      END
