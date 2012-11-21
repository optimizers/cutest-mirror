! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CDIMSJ( data, status, nnzj )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER, INTENT( OUT ) :: status
      INTEGER :: nnzj

!  Compute the space required to store the Jacobian matrix of the 
!  constraints/objective function of a problem initially written in 
!  Standard Input Format (SIF).

!  Based on the minimization subroutine data%laNCELOT/SBMIN
!  by Conn, Gould and Toint.

!  Nick Gould, for CGT productions,
!  August 1999.


! ---------------------------------------------------------------------




! ---------------------------------------------------------------------




! ---------------------------------------------------------------------


! ---------------------------------------------------------------------

!  integer variables from the GLOBAL common block.


!  integer variables from the LOCAL common block.




!  local variables.

      INTEGER :: ig

!  The total space is stored in nnzj

      nnzj = 0

!  Allow space for constraint groups

      DO 10 ig = 1, data%ng
         IF ( data%KNDOFC( ig ) /= 0 ) nnzj = nnzj + &
              data%IWORK( data%lstagv + ig + 1 ) - data%IWORK( data%lstagv + ig )
   10 CONTINUE

!  Add space for the (dense) gradient of the objective function.

      nnzj = nnzj + data%numvar
      status = 0
      RETURN

!  end of CDIMSJ.

      END
