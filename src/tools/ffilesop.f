      SUBROUTINE FORTRAN_OPEN( FUNIT, FNAME, IERR )

!     Open file FNAME using unit number FUNIT. 
!     If successful, IERR = 0 is returned, otherwise IERR = 1.
!     FORTRAN_OPEN() is particularly intended to be called from C
!     when a unit number is required.

      IMPLICIT NONE
      INTEGER :: FUNIT, IERR
      CHARACTER*64 FNAME

      OPEN( FUNIT, FILE=FNAME, STATUS='UNKNOWN', ERR=9000 )
      IERR = 0
      RETURN

 9000 IERR = 1
      RETURN
      END

      SUBROUTINE FORTRAN_CLOSE( FUNIT, IERR )

!     Close a stream unit previously opened by FORTRAN_OPEN
!     Exit value: 0 = successful return, 1 = error.

      IMPLICIT NONE
      INTEGER :: FUNIT, IERR

      CLOSE( FUNIT, ERR = 9001 )
      IERR = 0
      RETURN

 9001 IERR = 1
      RETURN
      END
