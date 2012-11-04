! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CONNAMES( data, m, GNAME )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: m
      CHARACTER ( LEN = 10 ) :: GNAME( m )

!  Obtain the names of the general constraints.

!  D. Orban, for GOT productions, Aug 2005, from the version of
!  Nick Gould, for CGT productions.
!  September 1992.


! ---------------------------------------------------------------------




! ---------------------------------------------------------------------



! ---------------------------------------------------------------------


! ---------------------------------------------------------------------

!  integer variables from the GLOBAL common block.


!  integer variables from the LOCAL common block.


!  local variables.

      INTEGER :: i, ig

!  Set the names of the general constraints.

!                            only if there are constraints.
      IF ( data%numcon > 0 ) THEN
         DO 20 ig = 1, data%ng
            IF ( data%KNDOFC( ig ) /= 0 ) &
                 GNAME( data%KNDOFC( ig ) ) = data%GNAMES( ig )
   20    CONTINUE
      END IF
      RETURN

!  end of CONNAMES.

      END
