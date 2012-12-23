! THIS VERSION: CUTEST 1.0 - 20/11/2012 AT 13:15 GMT.

!-*-*-*-*-*-  C U T E S T    C V A R T Y P E    S U B R O U T I N E  -*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 77 version originally released in CUTEr, December 1999
!   fortran 2003 version released in CUTEst, 20th November 2012

      SUBROUTINE CUTEST_cvartype( data, status, n, X_type )
      USE CUTEST

!  dummy arguments

      TYPE ( CUTEST_data_type ), INTENT( IN ) :: data
      INTEGER, INTENT( IN ) :: n
      INTEGER, INTENT( OUT ) :: status
      INTEGER, INTENT( OUT ), DIMENSION( n ) :: X_type

!  --------------------------------------------------------------
!  determine the type (continuous, 0-1, integer) of each variable
!  --------------------------------------------------------------

!  set the type of each variable (0 = continuous, 1 = 0-1, 2 = integer)

      X_type( : n ) = data%ITYPEV( : n )
      status = 0
      RETURN

!  end of subroutine CUTEST_cvartype

      END SUBROUTINE CUTEST_cvartype
