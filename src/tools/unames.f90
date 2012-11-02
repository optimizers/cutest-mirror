! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE UNAMES( data, N, PNAME, VNAME )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: N
      CHARACTER ( LEN = 10 ) :: PNAME, VNAME( N )

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

      INTEGER :: I

!  Set the problem name.

      PNAME = data%VNAMES( N + 1 )

!  Set the names of the variables.

      DO 10 I = 1, N
        VNAME( I ) = data%VNAMES( I )
   10 CONTINUE
      RETURN

!  end of UNAMES.

      END
