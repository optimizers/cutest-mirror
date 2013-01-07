!     ( Last modified on 3 Jan 2013 at 13:20:00 )
Program COBMA

!  COBYLA test driver for problems derived from SIF files.

!  A. R. Conn and Ph. Toint (based upon Nick Gould's vf13ma.f)
!  January 1995.
!  Revised for CUTEst, Nick Gould, January 2013

!  Fortran 90/95 version, D. Orban, December 2006

  Use CUTEst_problem

  Implicit None

  INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
  Type( CUTEST_problem_type ) :: prob

  Integer :: M, MAXFUN, LW, LIW, status
  Integer :: IPRINT, I, mgeq, NFIX
  Integer, Dimension(:), Allocatable :: IW

  Real( Kind = wp ), Dimension(:), Allocatable :: W
  Real( Kind = wp ) :: RHOBEG, RHOEND, F
  Real( Kind = wp ), Dimension(2) :: CPU
  Real( Kind = wp ), Dimension(7) :: CALLS
  Integer :: ierr
  INTEGER :: io_buffer = 11
  Integer, Parameter :: INPUT = 55, INDR = 46, IOUT = 6

!  Open the relevant file.

  Open( INPUT, FILE = 'OUTSDIF.d', FORM = 'FORMATTED', STATUS = 'OLD' )
  Rewind INPUT

!  Initialize problem data structure

  prob%allocate_H = .False.     ! No need for Hessian of objective/Lagrangian
  prob%allocate_J = .False.     ! No need for Jacobian of constraints
  Call CUTEST_problem_setup( status, prob, INPUT )
  IF ( status /= 0 ) GO TO 910

!  Set up the data structures necessary to hold the problem functions.

  Call CUTEST_csetup( status, INPUT, IOUT, io_buffer, prob%n, prob%m, prob%x,  &
                      prob%x_l, prob%x_u, prob%y, prob%c_l, prob%c_u,          &
                      prob%equation, prob%linear, .True., .False., .False. )
  IF ( status /= 0 ) GO TO 910
  Close( INPUT )

!  Allocate temporary work arrays

  lw  = prob%n * ( 3 * prob%n + 2 * prob%m + 11 ) + 4 * prob%m + 6
  Allocate( W( lw ), STAT = status )
  If ( status /= 0 ) GO TO 990
  liw = prob%n + 1
  Allocate( IW( liw ), STAT = status )
  If ( status /= 0 ) GO TO 990

!  Count the number of general equality constraints and ignore them by
!  shifting the remaining constraints at the beginning of the constraint list

  mgeq = 0
  Do I = 1, prob%m
     If ( prob%equation( I ) ) mgeq = mgeq + 1
  End Do
  If ( mgeq > 0 ) Then
     Write( 6, 3090 ) mgeq
     Do I = prob%m, mgeq + 1, -1
        prob%c_u( I - mgeq ) = prob%c_u( I )
        prob%c_l( I - mgeq ) = prob%c_l( I )
     End Do
  End If
  M = prob%m - mgeq

!  If constraints have both lower and upper bounds, they have to be
!  included twice!

  Do I = 1, prob%m - mgeq
     If ( prob%c_l( I ) > -INFTY .And. prob%c_u( I ) < INFTY ) M = M + 1
  End Do

!  Include any simple bounds.

  NFIX = 0
  Do I = 1, prob%n
     If ( prob%x_l( I ) ==  prob%x_u( I ) ) Then
        NFIX = NFIX + 1
     Else
        If ( prob%x_l( I ) > -INFTY ) M = M + 1
        If ( prob%x_u( I ) <  INFTY ) M = M + 1
     End If
  End Do
  If ( NFIX > 0 ) Write( 6, 3020 ) NFIX

!  Open the Spec file for the method.

  Open( INDR, FILE = 'COBYLA.SPC', FORM = 'FORMATTED', STATUS = 'OLD')
  Rewind INDR

!  Read input Spec data.

!  RHOBEG = the size of the simplex initially
!  RHOEND = the size of the simplex at termination
!  MAXFUN = the maximum number of function calls allowed.
!  IPRINT   should be set to 0, 1, 2 or 3, it controls the amount of printing

!  Set up algorithmic input data.

  Read ( INDR, 1000 ) RHOBEG, RHOEND, MAXFUN, IPRINT
  Close ( INDR )

!  Evaluate the objective function and constraints.

  Call CALCFC( prob, F, mgeq )

!  Perform the minimization.

  Call COBYLA( prob%n, M, prob%x, RHOBEG, RHOEND, IPRINT, MAXFUN, W, IW)

!  Output report

  Call CUTEST_creport( status, CALLS, CPU )
  IF ( status /= 0 ) GO TO 910

  Call CUTEST_cnames( status, prob%n, prob%m, prob%pname, prob%vnames, prob%cnames )
  Call CALCFC( prob, F, mgeq )
  Write( 6, 2110 ) ( I, prob%vnames( I ), prob%x( I ), &
       prob%x_l( I ), prob%x_u( I ), I = 1, prob%n )
  If ( prob%m > 0 ) Write( 6, 2120 ) ( I, prob%cnames( I ), prob%c( I ), &
       prob%c_l( I ), prob%c_u( I ), prob%linear( I ), I = 1, prob%m )
  Write( 6, 2000 ) prob%pname, prob%n, prob%m, CALLS(1), CALLS(5), F, &
       CPU(1), CPU(2)

!  Clean-up data structures

  Call CUTEST_problem_terminate( status, prob )
  IF ( status /= 0 ) GO TO 910
  Deallocate( IW, STAT = ierr )
  Deallocate( W, STAT = ierr )
  Stop

  910 CONTINUE
      WRITE( 6, "( ' CUTEst error, status = ', i0, ', stopping' )") status
      STOP

  990 CONTINUE
      WRITE( out, "( ' Allocation error, status = ', I0 )" ) status
      STOP

!  Non-executable statements

2000 Format( /, 24('*'), ' CUTEst statistics ', 24('*') //, &
          ' Code used               :  COBYLA ',  /,    &
          ' Problem                 :  ', A10,    /,    &
          ' # variables             =      ', I10 /,    &
          ' # constraints           =      ', I10 /,    &
          ' # objective functions   =        ', F8.2 /, &
          ' # constraints functions =        ', F8.2 /, &
          ' Final f                 = ', E15.7 /,       &
          ' Set up time             =      ', 0P, F10.2, ' seconds' /, &
          ' Solve time              =      ', 0P, F10.2, ' seconds' //, &
          66('*') / )
1000 Format( D12.4, /, D12.4, /,I6, /, I6 )
2110 Format( /, ' The variables:', /, &
          '     I name          value    lower bound upper bound', &
          /, ( I6, 1X, A10, 1P, 3D12.4 ) )
2120 Format( /, ' The constraints:', /, &
          '     I name          value    lower bound upper bound', &
          ' linear? ', &
          /, ( I6, 1X, A10, 1P, 3D12.4, 5X, L1 ) )
3000 Format( /,'  ** Program CSETUP: array length ', A6, ' too small.', &
          /,'  -- Miminimization abandoned.', &
          /,'  -- Increase the parameter ', A6, ' by at least ', I8, &
            ' and restart.'  )
3020 Format( /,'  ** Warning from COBMA. **', &
          /,'     In the problem as stated , ', I6, &
            ' variables are fixed: they are changed to free.' )
3090 Format( /,'  ** Warning from COBMA. **', &
          /,'     The problem as stated includes ', I6, &
            ' equality constraints: they are ignored ' )

!  End of COBMA.

End Program COBMA
Subroutine CALCFC( prob, F, mgeq )

!  Evaluates the objective function value in a format compatible with COBYLA,
!  but using the CUTEst tools.

!  A. R. Conn and Ph. Toint
!  January 1995.

!  Fortran 90/95 version, D. Orban, December 2006

  Use CUTEst_problem
  Implicit None

  INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
  Type( CUTEst_prob_type ), Intent( Inout ) :: prob
  Real( Kind = wp ), Intent( Inout ) :: F
  Integer, Intent( In ) :: mgeq
  Integer :: I, MT, status

!  Evaluate the objective function and constraints.

  Call CUTEST_cfn( status, prob%n, prob%m, prob%x, F, prob%c )
  IF ( status /= 0 ) GO TO 910

!  If there are equality constraints, ignore them
!  and shift all the inequality constraint values.

  If ( mgeq > 0 ) Then
     Do I = prob%m, mgeq + 1, - 1
        prob%c( I - mgeq ) = prob%c( I )
     End Do
  End If

!  If constraints have both lower and upper bounds, they have to
!  be included twice! Reverse the signs of less-than-or-equal-to constraints.

  MT = prob%m - mgeq + 1
  Do I = 1, prob%m - mgeq
     If ( prob%c_l( I ) > -INFTY .And. prob%c_u( I ) < INFTY ) Then
        prob%c( I )  = prob%c_u( I ) - prob%c( I )
        prob%c( MT ) = prob%c( I ) - prob%c_l( I )
        MT      = MT + 1
     Else If ( prob%c_l( I ) > -INFTY ) Then
        prob%c( I )  = prob%c( I ) - prob%c_l( I )
     Else If ( prob%c_u( I ) < INFTY ) Then
        prob%c( I )  = prob%c_u( I ) - prob%c( I )
     End If
  End Do

!  Include any simple bounds, including fixed variables.

  Do I = 1, prob%n
     If ( prob%x_l( I ) /=  prob%x_u( I ) ) Then
        If ( prob%x_l( I ) > -INFTY ) Then
           prob%c( MT ) = prob%x( I ) - prob%x_l( I )
           MT      = MT + 1
        End If
        If ( prob%x_u( I ) < INFTY ) Then
           prob%c( MT ) = prob%x_u( I ) - prob%x( I )
           MT      = MT + 1
        End If
     End If
  End Do
  Return

  910 CONTINUE
      WRITE( 6, "( ' CUTEst error, status = ', i0, ', stopping' )") status
      STOP

!  End of CALCFC.

End Subroutine CALCFC
