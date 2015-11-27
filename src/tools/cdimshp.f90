! THIS VERSION: CUTEST 1.3 - 24/11/2015 AT 15:05 GMT.

!-*-*-*-*-*-*-  C U T E S T    C D I M S H P    S U B R O U T I N E  -*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 2003 version released in CUTEst, 24th November 2015

      SUBROUTINE CUTEST_cdimshp( status, nnzshp )
      USE CUTEST

!  dummy arguments

      INTEGER, INTENT( OUT ) :: nnzshp, status

!  -----------------------------------------------------------------
!  compute the space required to store the matrix of products of the 
!  constraint Hessians with a vector of a problem initially written 
!  in Standard Input Format (SIF)
!  -----------------------------------------------------------------

      CALL CUTEST_cdimshp_threadsafe( CUTEST_data_global, status, nnzshp )
      RETURN

!  end of sunroutine CUTEST_cdimshp

      END SUBROUTINE CUTEST_cdimshp

!-  C U T E S T   C D I M S H P _ t h r e a d s a f e   S U B R O U T I N E  -

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 2003 version released in CUTEst, 24th November 2015

      SUBROUTINE CUTEST_cdimshp_threadsafe( data, status, nnzshp )
      USE CUTEST

!  dummy arguments

      TYPE ( CUTEST_data_type ), INTENT( IN ) :: data
      INTEGER, INTENT( OUT ) :: nnzshp, status

!  -----------------------------------------------------------------
!  compute the space required to store the matrix of products of the 
!  constraint Hessians with a vector of a problem initially written 
!  in Standard Input Format (SIF)
!  -----------------------------------------------------------------

!  local variables

      INTEGER :: ig

!  the total space is stored in nnzshp

      nnzshp = 0

!  allow space for constraint groups

      DO ig = 1, data%ng
        IF ( data%KNDOFC( ig ) /= 0 )                                         &
          nnzshp = nnzshp + data%ISTAGV( ig + 1 ) - data%ISTAGV( ig )
      END DO

      status = 0
      RETURN

!  end of sunroutine CUTEST_cdimshp_threadsafe

      END SUBROUTINE CUTEST_cdimshp_threadsafe
