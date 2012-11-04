! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE VARNAMES( data, n, VNAME )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n
      CHARACTER ( LEN = 10 ) :: VNAME( n )

!  Obtain the names of the problem variables.

!  D. Orban, for GOT productions, Aug 2005, from the version of
!  Nick Gould, for CGT productions.
!  September 1992.

!  local variables.

      INTEGER :: i

!  Set the names of the variables.

      DO 10 i = 1, n
        VNAME( i ) = data%VNAMES( i )
   10 CONTINUE
      RETURN

!  end of VARNAMES.

      END
