      SUBROUTINE FORTRAN_OPEN( funit, FNAME, ierr )

!     Open file FNAME using unit number FUNIT. 
!     If successful, ierr = 0 is returned, otherwise ierr = 1.
!     FORTRAN_OPEN() is particularly intended to be called from C
!     when a unit number is required.

      IMPLICIT NONE
      INTEGER :: funit, ierr
      CHARACTER*64 FNAME

      OPEN( funit, FILE=FNAME, STATUS='UNKNOWN', ERR=9000 )
      ierr = 0
      RETURN

 9000 ierr = 1
      RETURN
      END

      SUBROUTINE FORTRAN_CLOSE( funit, ierr )

!     Close a stream unit previously opened by FORTRAN_OPEN
!     Exit value: 0 = successful return, 1 = error.

      IMPLICIT NONE
      INTEGER :: funit, ierr

      CLOSE( funit, ERR = 9001 )
      ierr = 0
      RETURN

 9001 ierr = 1
      RETURN
      END
