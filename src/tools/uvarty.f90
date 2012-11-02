! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE UVARTY( data, N, IVARTY )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: N
      INTEGER :: IVARTY( N )

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

      INTEGER :: I

!  Set the type of each variable.

      DO 10 I = 1, N
        IVARTY( I ) = data%ITYPEV( I )
   10 CONTINUE
      RETURN

!  end of UVARTY.

      END
