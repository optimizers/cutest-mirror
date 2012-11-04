! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE UVARTY( data, n, IVARTY )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n
      INTEGER :: IVARTY( n )

!  Determine the type (continuous, 0-1, integer) of each variable

!  Nick Gould, for CGT productions.
!  December 1999.


! ---------------------------------------------------------------------




! ---------------------------------------------------------------------



! ---------------------------------------------------------------------


! ---------------------------------------------------------------------

!  integer variables from the GLOBAL common block.


!  integer variables from the LOCAL common block.


!  local variables.

      INTEGER :: i

!  Set the type of each variable.

      DO 10 i = 1, n
        IVARTY( i ) = data%ITYPEV( i )
   10 CONTINUE
      RETURN

!  end of UVARTY.

      END
