! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CNAMES( data, N, M, PNAME, VNAME, GNAME )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: N, M
      CHARACTER ( LEN = 10 ) :: PNAME, VNAME( N ), GNAME( M )

!  Obtain the names of the problem, its variables and
!  general constraints.

!  Nick Gould, for CGT productions.
!  September 1992.


! ---------------------------------------------------------------------




! ---------------------------------------------------------------------



! ---------------------------------------------------------------------


! ---------------------------------------------------------------------

!  integer variables from the GLOBAL common block.


!  integer variables from the LOCAL common block.


!  local variables.

      INTEGER :: I, IG

!  Set the problem name.

      PNAME = data%VNAMES( N + 1 )

!  Set the names of the variables.

      DO 10 I = 1, N
        VNAME( I ) = data%VNAMES( I )
   10 CONTINUE

!  Set the names of the general constraints.

!                            only if there are constraints.
      IF ( data%numcon > 0 ) THEN
         DO 20 IG = 1, data%ng
            IF ( data%KNDOFC( IG ) /= 0 ) &
                 GNAME( data%KNDOFC( IG ) ) = data%GNAMES( IG )
   20    CONTINUE
      END IF
      RETURN

!  end of CNAMES.

      END
