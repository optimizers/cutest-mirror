! THIS VERSION: CUTEST 1.0 - 04/11/2012 AT 12:50 GMT.

!-*-*-*-*-*-*-  C U T E S T    C R E P R T    S U B R O U T I N E  -*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal authors: Nick Gould and Philippe Toint

!  History -
!   fortran 77 version originally released in CUTEr, 23rd December, 2000
!   fortran 2003 version released in CUTEst, 4th November 2012

      SUBROUTINE CREPRT( data, status, CALLS, TIME )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  Dummy arguments

      INTEGER, INTENT( OUT ) :: status
      REAL ( KIND = wp ), DIMENSION( 7 ):: CALLS
      REAL ( KIND = wp ), DIMENSION( 2 ):: TIME

!  ------------------------------------------------------------------------
!  return the values of counters maintained by the CUTEst tools. 
!  The counters are:

!    CALLS( 1 ): number of calls to the objective function
!    CALLS( 2 ): number of calls to the objective gradient
!    CALLS( 3 ): number of calls to the objective Hessian
!    CALLS( 4 ): number of Hessian times vector products
!    CALLS( 5 ): number of calls to the constraint functions
!    CALLS( 6 ): number of calls to the constraint gradients
!    CALLS( 7 ): number of calls to the constraint Hessians

!    TIME( 1 ): CPU time (in seconds) for CSETUP
!    TIME( 2 ): CPU time (in seconds) since the end of CSETUP

!  Note that each constraint function is counted separately. Evaluating all
!  the constraints thus results in data%pnc evaluations, where data%pnc is 
!  the number of constraints in the problem.  Note that data%pnc does not
!  include repetitions for constraints having full ranges
!  ------------------------------------------------------------------------

!  local variable

      REAL ( KIND = wp ) :: time_now

      CALL CPU_TIME( time_now )

      TIME( 1 ) = data%sutime
      TIME( 2 ) = time_now - data%sttime

      CALLS( 1 ) = data%nc2of
      CALLS( 2 ) = data%nc2og
      CALLS( 3 ) = data%nc2oh
      CALLS( 4 ) = data%nhvpr
      IF( data%pnc > 0 ) THEN
        CALLS( 5 ) = data%nc2cf / data%pnc
        CALLS( 6 ) = data%nc2cg / data%pnc
        CALLS( 7 ) = data%nc2ch / data%pnc
      ELSE
        CALLS( 5 ) = data%nc2cf
        CALLS( 6 ) = data%nc2cg
        CALLS( 7 ) = data%nc2ch
      ENDIF

      status = 0
      RETURN

!  End of subroutine CREPRT

      END SUBROUTINE CREPRT
