!   ( Last modified on 3 Jan 2013 at 13:20:00 )

      PROGRAM COBYLA_main

!  COBYLA test driver for problems derived from SIF files

!  A. R. Conn and Ph. Toint (based upon Nick Gould's vf13ma.f)
!  January 1995.
!  Fortran 90/95 version, D. Orban, December 2006
!  Revised for CUTEst, Nick Gould, January 2013

      USE CUTEst_problem

      IMPLICIT NONE
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      TYPE( CUTEST_problem_type ) :: prob
      INTEGER :: m, maxfun, lw, liw, status, iprint, i, mgeq, nfix, ierr
      INTEGER, DIMENSION(:), ALLOCATABLE :: IW
      REAL( KIND = wp ) :: rhobeg, rhoend, f
      REAL( KIND = wp ), PARAMETER :: infty = 1.0D+19
      REAL( KIND = wp ), DIMENSION(:), ALLOCATABLE :: W
      REAL( KIND = wp ), DIMENSION( 2 ) :: CPU
      REAL( KIND = wp ), DIMENSION( 7 ) :: CALLS
      INTEGER :: io_buffer = 11
      INTEGER, PARAMETER :: input = 55, indr = 46, out = 6

!  open the relevant file

      OPEN( input, FILE = 'OUTSDIF.d', FORM = 'FORMATTED', STATUS = 'OLD' )
      REWIND( input )

!  initialize problem data structure

      prob%allocate_H = .FALSE.   ! No need for Hessian of objective/Lagrangian
      prob%allocate_J = .FALSE.   ! No need for Jacobian of constraints
      CALL CUTEST_problem_setup( status, prob, INPUT )
      IF ( status /= 0 ) GO TO 910

!  set up the data structures necessary to hold the problem functions.

      CALL CUTEST_csetup( status, input, out, io_buffer, prob%n, prob%m,       &
                          prob%x, prob%x_l, prob%x_u, prob%y, prob%c_l,        &
                          prob%c_u, prob%equation, prob%linear,                &
                          .TRUE., .FALSE., .FALSE. )
      IF ( status /= 0 ) GO TO 910
      CLOSE( input )

!  allocate temporary work arrays

      lw  = prob%n * ( 3 * prob%n + 2 * prob%m + 11 ) + 4 * prob%m + 6
      ALLOCATE( W( lw ), STAT = status )
      IF ( status /= 0 ) GO TO 990
      liw = prob%n + 1
      ALLOCATE( IW( liw ), STAT = status )
      IF ( status /= 0 ) GO TO 990

!  count the number of general equality constraints and ignore them by
!  shifting the remaining constraints at the beginning of the constraint list

      mgeq = 0
      DO i = 1, prob%m
        IF ( prob%equation( i ) ) mgeq = mgeq + 1
      END DO
      IF ( mgeq > 0 ) THEN
        WRITE( 6, 3090 ) mgeq
        DO i = prob%m, mgeq + 1, -1
          prob%c_u( i - mgeq ) = prob%c_u( i )
          prob%c_l( i - mgeq ) = prob%c_l( i )
        END DO
      END IF
      m = prob%m - mgeq

!  if constraints have both lower and upper bounds, they must be included twice!

      DO i = 1, prob%m - mgeq
         IF ( prob%c_l( i ) > - infty .And. prob%c_u( i ) < infty ) m = m + 1
      END DO

!  include any simple bounds

      nfix = 0
      DO i = 1, prob%n
        IF ( prob%x_l( i ) ==  prob%x_u( i ) ) THEN
          nfix = nfix + 1
        Else
          IF ( prob%x_l( i ) > -infty ) m = m + 1
          IF ( prob%x_u( i ) <  infty ) m = m + 1
        END IF
      END DO
      IF ( nfix > 0 ) WRITE( 6, 3020 ) NFIX

!  open the Spec file for the method

      OPEN( indr, FILE = 'COBYLA.SPC', FORM = 'FORMATTED', STATUS = 'OLD')
      REWIND( indr )

!  read input Spec data

!  RHOBEG = the size of the simplex initially
!  RHOEND = the size of the simplex at termination
!  MAXFUN = the maximum number of function calls allowed.
!  IPRINT   should be set to 0, 1, 2 or 3, it controls the amount of printing

!  set up algorithmic input data

      READ ( indr, 1000 ) RHOBEG, RHOEND, MAXFUN, IPRINT
      CLOSE ( indr )

!  evaluate the objective function and constraints

      CALL COBYLA_evalf( prob, f, mgeq )

!  perform the minimization

      CALL COBYLA( prob%n, m, prob%x, rhobeg, rhoend, iprint, maxfun, W, IW )

!  output report

      CALL CUTEST_creport( status, CALLS, CPU )
      IF ( status /= 0 ) GO TO 910

      CALL CUTEST_cnames( status, prob%n, prob%m, prob%pname,                  &
                          prob%vnames, prob%cnames )
      CALL COBYLA_evalf( prob, f, mgeq )
      WRITE( out, 2110 ) ( i, prob%vnames( i ), prob%x( i ),                   &
           prob%x_l( i ), prob%x_u( i ), i = 1, prob%n )
      IF ( prob%m > 0 ) WRITE( 6, 2120 ) ( I, prob%cnames( i ), prob%c( i ),   &
           prob%c_l( i ), prob%c_u( i ), prob%linear( i ), i = 1, prob%m )
      WRITE( out, 2000 ) prob%pname, prob%n, prob%m, CALLS( 1 ), CALLS( 5 ),   &
           f, CPU( 1 ), CPU( 2 )

!  clean-up data structures

      CALL CUTEST_problem_terminate( status, prob )
      IF ( status /= 0 ) GO TO 910
      DEALLOCATE( IW, STAT = ierr )
      DEALLOCATE( W, STAT = ierr )
      CALL CUTEST_cterminate( status )
      STOP

!  error returns

  910 CONTINUE
      WRITE( 6, "( ' CUTEst error, status = ', i0, ', stopping' )") status
      STOP

  990 CONTINUE
      WRITE( out, "( ' Allocation error, status = ', I0 )" ) status
      STOP

!  Non-executable statements

2000 FORMAT( /, 24('*'), ' CUTEst statistics ', 24('*') //,                    &
          ' Code used               :  COBYLA ',  /,                           &
          ' Problem                 :  ', A10,    /,                           &
          ' # variables             =      ', I10 /,                           &
          ' # constraints           =      ', I10 /,                           &
          ' # objective functions   =        ', F8.2 /,                        &
          ' # constraints functions =        ', F8.2 /,                        &
          ' Final f                 = ', E15.7 /,                              &
          ' Set up time             =      ', 0P, F10.2, ' seconds' /,         &
          ' Solve time              =      ', 0P, F10.2, ' seconds' //,        &
          66('*') / )
1000 FORMAT( D12.4, /, D12.4, /,I6, /, I6 )
2110 FORMAT( /, ' The variables:', /, &
          '     i name          value    lower bound upper bound',             &
          /, ( I6, 1X, A10, 1P, 3D12.4 ) )
2120 FORMAT( /, ' The constraints:', /, &
          '     i name          value    lower bound upper bound',             &
          ' linear? ', &
          /, ( I6, 1X, A10, 1P, 3D12.4, 5X, L1 ) )
3000 FORMAT( /,'  ** Program CSETUP: array length ', A6, ' too small.',        &
          /,'  -- Miminimization abandoned.', &
          /,'  -- Increase the parameter ', A6, ' by at least ', I0,           &
            ' and restart.'  )
3020 FORMAT( /,'  ** Warning from COBYLA_main. **',                            &
          /,'     In the problem as stated , ', I0,                            &
            ' variables are fixed: they are changed to free.' )
3090 FORMAT( /,'  ** Warning from COBYLA_main. **',                            &
          /,'     The problem as stated includes ', I0,                        &
            ' equality constraints: they are ignored ' )

!  End of COBYLA_main

      END PROGRAM COBYLA_main

      SUBROUTINE COBYLA_evalf( prob, f, mgeq )

!  evaluates the objective function value in a format compatible with COBYLA
!  using the CUTEst tools.

     USE CUTEst_problem
     IMPLICIT NONE

     INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
     REAL( KIND = wp ), PARAMETER :: infty = 1.0D+19
     TYPE( CUTEst_problem_type ), INTENT( INOUT ) :: prob
     REAL( Kind = wp ), INTENT( INOUT ) :: f
     INTEGER, INTENT( IN ) :: mgeq
     INTEGER :: i, mt, status

!  Evaluate the objective function and constraints.

      CALL CUTEST_cfn( status, prob%n, prob%m, prob%x, F, prob%c )
      IF ( status /= 0 ) GO TO 910

!  IF there are equality constraints, ignore them
!  and shift all the inequality constraint values.

      IF ( mgeq > 0 ) THEN
        DO i = prob%m, mgeq + 1, - 1
          prob%c( i - mgeq ) = prob%c( i )
        END DO
      END IF

!  IF constraints have both lower and upper bounds, they have to
!  be included twice! Reverse the signs of less-than-or-equal-to constraints.

      mt = prob%m - mgeq + 1
      DO i = 1, prob%m - mgeq
        IF ( prob%c_l( i ) > - infty .And. prob%c_u( i ) < infty ) THEN
          prob%c( i )  = prob%c_u( i ) - prob%c( i )
          prob%c( MT ) = prob%c( i ) - prob%c_l( i )
          mt = mt + 1
        ELSE IF ( prob%c_l( i ) > - infty ) THEN
          prob%c( i )  = prob%c( i ) - prob%c_l( i )
        ELSE IF ( prob%c_u( i ) < infty ) THEN
          prob%c( i )  = prob%c_u( i ) - prob%c( i )
        END IF
      END DO

!  Include any simple bounds, including fixed variables.

      DO i = 1, prob%n
        IF ( prob%x_l( i ) /=  prob%x_u( i ) ) THEN
          IF ( prob%x_l( i ) > - infty ) THEN
            prob%c( mt ) = prob%x( i ) - prob%x_l( i )
            mt = mt + 1
          END IF
          IF ( prob%x_u( i ) < infty ) THEN
            prob%c( mt ) = prob%x_u( i ) - prob%x( i )
            mt = mt + 1
          END IF
        END IF
      END DO
      RETURN

  910 CONTINUE
      WRITE( 6, "( ' CUTEst error, status = ', i0, ', stopping' )") status
      STOP

!  End of COBYLA_evalf.

      END SUBROUTINE COBYLA_evalf
