! THIS VERSION: CUTEST 1.0 - 19/11/2012 AT 13:10 GMT.

!-*-*-*-*-*-  C U T E S T    C O N N A M E S    S U B R O U T I N E  -*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Dominique Orban

!  History -
!   fortran 77 version originally released in CUTEr, August 2005
!   fortran 2003 version released in CUTEst, 19th November 2012

      SUBROUTINE CONNAMES( data, status, m, CNAME )
      USE CUTEST

!  dummy arguments

      TYPE ( CUTEST_data_type ), INTENT( INOUT ) :: data
      INTEGER, INTENT( IN ) :: m
      INTEGER, INTENT( OUT ) :: status
      CHARACTER ( LEN = 10 ), INTENT( OUT ), DIMENSION( m ) :: CNAME

!  -------------------------------------------
!  obtain the names of the general constraints
!  -------------------------------------------

      INTEGER :: ig

!  Set the names of the general constraints

      IF ( data%numcon > 0 ) THEN
        DO ig = 1, data%ng
          IF ( data%KNDOFC( ig ) /= 0 )                                        &
             CNAME( data%KNDOFC( ig ) ) = data%GNAMES( ig )
        END DO
      END IF
      status = 0
      RETURN

!  end of subroutine CONNAMES

      END SUBROUTINE CONNAMES
