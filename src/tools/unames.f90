! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE UNAMES( data, n, PNAME, VNAME )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n
      CHARACTER ( LEN = 10 ) :: PNAME, VNAME( n )

!  Obtain the names of the problem and its variables.

!  Nick Gould, for CGT productions.
!  September 1992.


! ---------------------------------------------------------------------




! ---------------------------------------------------------------------



! ---------------------------------------------------------------------


! ---------------------------------------------------------------------

!  integer variables from the GLOBAL common block.


!  integer variables from the LOCAL common block.


!  local variables.

      INTEGER :: i

!  Set the problem name.

      PNAME = data%VNAMES( n + 1 )

!  Set the names of the variables.

      DO 10 i = 1, n
        VNAME( i ) = data%VNAMES( i )
   10 CONTINUE
      RETURN

!  end of UNAMES.

      END
