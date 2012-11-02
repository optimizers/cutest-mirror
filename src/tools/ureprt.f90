! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE UREPRT( data, CALLS, TIME )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  This routine returns the value of the various counters maintained by the
!  CUTEr tools to the user.  The counters are:

!    CALLS( 1 ): number of calls to the objective function
!    CALLS( 2 ): number of calls to the objective gradient
!    CALLS( 3 ): number of calls to the objective Hessian
!    CALLS( 4 ): number of Hessian times vector products

!    TIME( 1 ): CPU time (in seconds) for USETUP
!    TIME( 2 ): CPU time ( in seconds) since the end of USETUP

!  Note that each constraint function is counted separately.
!  Evaluating all the constraints thus results in data%pnc evaluations, where
!  data%pnc is the number of constraints in the problem.  Note that data%pnc does not
!  include repetitions for constraints having full ranges.  Also note that

!  N. Gould, D. Orban and Ph. Toint for CUTEr, 2001.

!  output arguments

      REAL :: CALLS( 4 )
      REAL :: TIME( 2 )



!  local variables

      REAL :: CPUTIM, DUM
      EXTERNAL :: CPUTIM

      TIME( 2 ) = CPUTIM( DUM ) - data%sttime
      TIME( 1 ) = data%sutime

      CALLS( 1 ) = data%nc2of
      CALLS( 2 ) = data%nc2og
      CALLS( 3 ) = data%nc2oh
      CALLS( 4 ) = data%nhvpr

      RETURN
      END
