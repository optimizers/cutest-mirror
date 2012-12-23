! THIS VERSION: CUTEST 1.0 - 04/11/2012 AT 12:40 GMT.

!-*-*-*-*-*-*-  C U T E S T    U R E P O R T    S U B R O U T I N E  -*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal authors: Nick Gould and Philippe Toint

!  History -
!   fortran 77 version originally released as UREPRT in CUTEr, December, 2000
!   fortran 2003 version released in CUTEst, 4th November 2012

      SUBROUTINE CUTEST_ureport( data, work, status, CALLS, TIME )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      TYPE ( CUTEST_work_type ) :: work
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  Dummy arguments

      INTEGER, INTENT( OUT ) :: status
      REAL ( KIND = wp ), DIMENSION( 4 ):: CALLS
      REAL ( KIND = wp ), DIMENSION( 2 ):: TIME

!  ------------------------------------------------------------------------
!  return the values of counters maintained by the CUTEst tools. 
!  The counters are:

!    CALLS( 1 ): number of calls to the objective function
!    CALLS( 2 ): number of calls to the objective gradient
!    CALLS( 3 ): number of calls to the objective Hessian
!    CALLS( 4 ): number of Hessian times vector products

!    TIME( 1 ): CPU time (in seconds) for USETUP
!    TIME( 2 ): CPU time ( in seconds) since the end of USETUP
!  ------------------------------------------------------------------------

!  local variable

      REAL ( KIND = wp ) :: time_now

      CALL CPU_TIME( time_now )

      TIME( 1 ) = data%sutime
      TIME( 2 ) = time_now - data%sttime

      CALLS( 1 ) = work%nc2of
      CALLS( 2 ) = work%nc2og
      CALLS( 3 ) = work%nc2oh
      CALLS( 4 ) = work%nhvpr

      status = 0
      RETURN

!  End of subroutine CUTEST_ureport

      END SUBROUTINE CUTEST_ureport
