! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CREPRT( CALLS, TIME )

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
!  Evaluating all the constraints thus results in PNC evaluations, where
!  PNC is the number of constraints in the problem.  Note that PNC does not
!  include repetitions for constraints having full ranges.  Also note that

!  N. Gould, D. Orban & Ph. Toint for CUTEr, 2001.

      IMPLICIT NONE
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  output arguments

      REAL :: CALLS( 7 )
      REAL :: TIME( 2 )

!  variables from the PRFCTS common block

      INTEGER :: NC2OF, NC2OG, NC2OH, NC2CF, NC2CG, NC2CH
      INTEGER :: NHVPR, PNC
      REAL :: SUTIME, STTIME
      COMMON / PRFCTS /  NC2OF, NC2OG, NC2OH, NC2CF, NC2CG, NC2CH, &
                         NHVPR, PNC, SUTIME, STTIME
      SAVE             / PRFCTS /

!  local variables

      REAL :: CPUTIM, DUM
      EXTERNAL :: CPUTIM

      TIME( 2 ) = CPUTIM( DUM ) - STTIME
      TIME( 1 ) = SUTIME

      CALLS( 1 ) = NC2OF
      CALLS( 2 ) = NC2OG
      CALLS( 3 ) = NC2OH
      CALLS( 4 ) = NHVPR
      IF( PNC > 0 ) THEN
        CALLS( 5 ) = NC2CF / PNC
        CALLS( 6 ) = NC2CG / PNC
        CALLS( 7 ) = NC2CH / PNC
      ELSE
        CALLS( 5 ) = NC2CF
        CALLS( 6 ) = NC2CG
        CALLS( 7 ) = NC2CH
      ENDIF

      RETURN
      END
