! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CREPRT( data, CALLS, TIME )
      USE CUTEST
      IMPLICIT NONE
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      REAL :: CALLS( 7 )
      REAL :: TIME( 2 )

!  This routine returns the value of the various counters maintained by the
!  CUTEr tools to the user.  The counters are:

!    CALLS( 1 ): number of calls to the objective function
!    CALLS( 2 ): number of calls to the objective gradient
!    CALLS( 3 ): number of calls to the objective Hessian
!    CALLS( 4 ): number of Hessian times vector products
!    CALLS( 5 ): number of calls to the constraint functions
!    CALLS( 6 ): number of calls to the constraint gradients
!    CALLS( 7 ): number of calls to the constraint Hessians

!    TIME( 1 ): CPU time (in seconds) for CSETUP
!    TIME( 2 ): CPU time (in seconds) since the end of CSETUP

!  Note that each constraint function is counted separately.
!  Evaluating all the constraints thus results in data%pnc evaluations, where
!  data%pnc is the number of constraints in the problem.  Note that data%pnc does not
!  include repetitions for constraints having full ranges.  Also note that

!  N. Gould, D. Orban & Ph. Toint for CUTEr, 2001.

!  output arguments



!  local variables

      REAL :: CPUTIM, DUM
      EXTERNAL :: CPUTIM

      TIME( 2 ) = CPUTIM( DUM ) - data%sttime
      TIME( 1 ) = data%sutime

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

      RETURN
      END
