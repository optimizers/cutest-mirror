! THIS VERSION: CUTEST 1.0 - 19/11/2012 AT 13:25 GMT.

!-*-*-*-*-*-  C U T E S T    V A R N A M E S    S U B R O U T I N E  -*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Dominique Orban

!  History -
!   fortran 77 version originally released in CUTEr, August 2005
!   fortran 2003 version released in CUTEst, 19th November 2012

      SUBROUTINE VARNAMES( data, status, n, VNAME )
      USE CUTEST

!  dummy arguments

      TYPE ( CUTEST_data_type ), INTENT( INOUT ) :: data
      INTEGER, INTENT( IN ) :: n
      INTEGER, INTENT( OUT ) :: status
      CHARACTER ( LEN = 10 ), INTENT( OUT ), DIMENSION( n ) :: VNAME

!  -----------------------------------------
!  obtain the names of the problem variables
!  -----------------------------------------

!  set the names of the variables

      VNAME( : n ) = data%VNAMES( : n )
      status = 0
      RETURN

!  end of subroutine VARNAMES

      END SUBROUTINE VARNAMES
