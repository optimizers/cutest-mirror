! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE UDIMEN( INPUT, N )
      INTEGER :: INPUT, N

!  Compute the basic array dimension for the unconstrained optimization tools.

!  Nick Gould, for CGT productions,
!  26th August, 1999.

      REWIND INPUT
      READ( INPUT, 1001 ) N
      REWIND INPUT
      RETURN

!  Non-executable statements.

 1001 FORMAT( I8 )

!  End of UDIMEN.

      END
