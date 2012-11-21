! THIS VERSION: CUTEST 1.0 - 19/11/2012 AT 13:00 GMT.

!-*-*-*-*-*-*-  C U T E S T    P B N A M E     S U B R O U T I N E  -*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Dominique Orban

!  History -
!   fortran 77 version originally released in CUTEr, August 2005
!   fortran 2003 version released in CUTEst, 19th November 2012

      SUBROUTINE PBNAME( data, status, pname )
      USE CUTEST

!  Dummy arguments

      TYPE ( CUTEST_data_type ), INTENT( INOUT ) :: data
      INTEGER, INTENT( OUT ) :: status
      CHARACTER ( LEN = 10 ), INTENT( OUT ) :: pname

!  ------------------------------
!  Obtain the name of the problem
!  ------------------------------

      pname = data%pname
      status = 0
      RETURN

!  end of subroutine PBNAME

      END SUBROUTINE PBNAME
