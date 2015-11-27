!   ( Last modified on 6 Oct 2015 at 15:50:00 )

      PROGRAM RAL_NLLS_main
      USE ISO_C_BINDING
      USE NLLS_MODULE

!  RAL_NLLS test driver for problems derived from SIF files

!  Nick Gould, October 2015

      IMPLICIT NONE

      type, extends( params_base_type ) :: user_type
         ! still empty
      end type user_type

      INTEGER :: status, i, m, n
      INTEGER( c_int ) :: len_work_integer, len_work_real
      REAL( c_double ), PARAMETER :: infty = 1.0D+19
      REAL( c_double ), DIMENSION( : ), ALLOCATABLE :: X, X_l, X_u
      REAL( c_double ), DIMENSION( : ), ALLOCATABLE :: Y, C_l, C_u, F
      REAL( c_double ), DIMENSION( : ), ALLOCATABLE ::  Work_real
      INTEGER( c_int ), DIMENSION( : ), ALLOCATABLE ::  Work_integer
      type( user_type ), target :: params
      TYPE( NLLS_inform_type ) :: inform
      TYPE( NLLS_control_type ) :: control
      LOGICAL, DIMENSION( : ), ALLOCATABLE  :: EQUATN, LINEAR
      CHARACTER ( LEN = 10 ) :: pname
      CHARACTER ( LEN = 20 ) :: summary_file = REPEAT( ' ', 20 )
      CHARACTER ( LEN = 10 ), ALLOCATABLE, DIMENSION( : )  :: VNAMES, CNAMES
      REAL( c_double ), DIMENSION( 2 ) :: CPU
      REAL( c_double ), DIMENSION( 7 ) :: CALLS
      INTEGER :: io_buffer = 11
      INTEGER :: summary_unit, iores
      INTEGER, PARAMETER :: input = 55, indr = 46, out = 6
      LOGICAL :: filexx

!  Interface blocks

     INTERFACE
       SUBROUTINE eval_F( status, n, m, X, F, params )
         USE ISO_C_BINDING
         import :: params_base_type
         INTEGER ( c_int ), INTENT( OUT ) :: status
         INTEGER ( c_int ), INTENT( IN ) :: n, m
         REAL ( c_double ), DIMENSION( n ), INTENT( IN ) :: X
         REAL ( c_double ), DIMENSION( m ), INTENT( OUT ) :: F
         class( params_base_type ), intent(in) :: params
       END SUBROUTINE eval_F
     END INTERFACE

     INTERFACE
       SUBROUTINE eval_J( status, n, m, X, J, params )
         USE ISO_C_BINDING
         import :: params_base_type
         INTEGER ( c_int ), INTENT( OUT ) :: status
         INTEGER ( c_int ), INTENT( IN ) :: n, m
         REAL ( c_double ), DIMENSION( n ), INTENT( IN ) :: X
         REAL ( c_double ), DIMENSION( m * n ), INTENT( OUT ) :: J
         class( params_base_type ), intent(in) :: params
       END SUBROUTINE eval_J
     END INTERFACE

     INTERFACE
       SUBROUTINE eval_HF( status, n, m, X, F, H, params )
         USE ISO_C_BINDING
         import :: params_base_type
         INTEGER ( c_int ), INTENT( OUT ) :: status
         INTEGER ( c_int ), INTENT( IN ) :: n, m
         REAL ( c_double ), DIMENSION( n ), INTENT( IN ) :: X
         REAL ( c_double ), DIMENSION( m ), INTENT( IN ) :: F
         REAL ( c_double ), DIMENSION( n**n ), INTENT( OUT ) :: H
         class( params_base_type ), intent(in) :: params
       END SUBROUTINE eval_HF
     END INTERFACE

!  open the relevant file

      OPEN( input, FILE = 'OUTSDIF.d', FORM = 'FORMATTED', STATUS = 'OLD' )
      REWIND( input )

!  compute problem dimensions

      CALL CUTEST_cdimen( status, input, n, m )
      IF ( status /= 0 ) GO TO 910

!  allocate space 

      ALLOCATE( X( n ), X_l( n ), X_u( n ), Y( m ), C_l( m ), C_u( m ),        &
                EQUATN( m ), LINEAR( m ), STAT = status )
      IF ( status /= 0 ) GO TO 990

!  initialize problem data structure

!  set up the data structures necessary to hold the problem functions.

      CALL CUTEST_csetup( status, input, out, io_buffer, n, m,                 &
                          X, X_l, X_u, Y, C_l, C_u, EQUATN, LINEAR, 0, 0, 0 )
      IF ( status /= 0 ) GO TO 910
      CLOSE( input )

!  allocate more space 

      DEALLOCATE( X_l, X_u, Y, C_l, C_u, EQUATN, LINEAR )
      len_work_integer = 0
      len_work_real = m + n * ( m + n )
      ALLOCATE( Work_integer( len_work_integer ), Work_real( len_work_real ),  &
                STAT = status )
      IF ( status /= 0 ) GO TO 990

!  open the Spec file for the method

      OPEN( indr, FILE = 'RAL_NLLS.SPC', FORM = 'FORMATTED', STATUS = 'OLD')
      REWIND( indr )

!  read input Spec data

!  error = unit for error messages
!  out = unit for information
!  print_level = level of output desired (<=0 = nothing)
!  nlls_method = method used (1=dogleg, 2=AINT, 3=More-Sorensen)
!  model = model used (1=first order, 2=Newton)
!  initial_radius = initial TR radius
!  stop_g_absolute = absolute stopping tolerance
!  stop_g_relative = relative stopping tolerance
!  summary_unit = write a one line summary to this unit (-ve = don't write)
!  summary_file = file name for summary (20 chars max)

!  set up algorithmic input data

      READ ( indr, "( I6, 4( /, I6 ), 3( /, E12.0 ), /, I6, /, A )" )          &
        control%error, control%out, control%print_level,                       &
        control%nlls_method, control%model, control%initial_radius,            &
        control%stop_g_absolute, control%stop_g_relative,                      &
        summary_unit, summary_file
      CLOSE ( indr )

write(6,*) summary_unit, summary_file

      IF ( summary_unit > 0 ) THEN
        INQUIRE( FILE = summary_file, EXIST = filexx )
        IF ( filexx ) THEN
           OPEN( summary_unit, FILE = summary_file, FORM = 'FORMATTED',        &
               STATUS = 'OLD', IOSTAT = iores )
        ELSE
           OPEN( summary_unit, FILE = summary_file, FORM = 'FORMATTED',        &
                STATUS = 'NEW', IOSTAT = iores )
        END IF
        IF ( iores /= 0 ) THEN 
          write( out, "( ' IOSTAT = ', I0, ' when opening file ', A,           &
        &  '. Stopping ' )" ) iores, summary_file
          STOP
        END IF
        CALL CUTEST_probname( status, pname )
        WRITE( summary_unit, "( A10 )" ) pname
      END IF

!  call the minimizer

      CALL RAL_NLLS( n, m, X, eval_F, eval_J, eval_HF,                         &
                     params, inform, control )

      WRITE( out , "( A, I0, A, I0)") 'status = ', inform%status,              &
          '       iter = ', inform%iter
      IF ( status /= 0 ) GO TO 910

!  output report

      CALL CUTEST_creport( status, CALLS, CPU )
      IF ( status /= 0 ) GO TO 910

      ALLOCATE( F( m ), VNAMES( n ), CNAMES( m ), STAT = status )
      CALL CUTEST_cnames( status, n, m, pname, VNAMES, CNAMES )
      CALL eval_F( status, n, m, X, F, params)

      WRITE( out, 2110 ) ( i, VNAMES( i ), X( i ), i = 1, n )
!     WRITE( out, 2120 ) ( i, CNAMES( i ), F( i ), i = 1, m )
      WRITE( out, 2000 ) pname, n, m, inform%obj, INT( CALLS( 5 ) ),           &
        INT( CALLS( 6 ) ), INT( CALLS( 7 ) ), CPU( 1 ), CPU( 2 )

!  write summary if required

      IF ( summary_unit > 0 ) THEN
        BACKSPACE( summary_unit )
        WRITE( summary_unit, "( A10, 4I6, ES12.4 )" )                          &
          pname, n, m, inform%status, inform%iter, inform%obj
        CLOSE(  summary_unit )
      END IF

!  clean-up data structures

      DEALLOCATE( X, F, VNAMES, CNAMES, Work_integer, Work_real,               &
                  STAT = status )
      IF ( status /= 0 ) GO TO 910
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
          ' Package used            :  RAL_NLLS ',  /,                         &
          ' Problem                 :  ', A10,    /,                           &
          ' # variables             =  ', I0, /,                               &
          ' # residuals             =  ', I0, /,                               &
          ' Final f                 =', ES15.7 /,                              &
          ' # residual evaluations  =  ', I0, /,                               &
          ' # Jacobian evaluations  =  ', I0, /,                               &
          ' # Hessian evaluations   =  ', I0, /,                               &
          ' Set up time             =  ', 0P, F0.2, ' seconds' /,              &
          ' Solve time              =  ', 0P, F0.2, ' seconds' //,             &
           66('*') / )
2110 FORMAT( /, ' The variables:', /, &
          '     i name          value',  /, ( I6, 1X, A10, 1P, D12.4 ) )
!2120 FORMAT( /, ' The constraints:', /, '     i name          value',         &
!         /, ( I6, 1X, A10, 1P, D12.4 ) )

!  End of RAL_NLLS_main

      END PROGRAM RAL_NLLS_main

      SUBROUTINE eval_F( status, n, m, X, F, params )
      USE ISO_C_BINDING
      use :: nlls_module, only : params_base_type
      
      INTEGER ( c_int ), INTENT( OUT ) :: status
      INTEGER ( c_int ), INTENT( IN ) :: n, m
      REAL ( c_double ), DIMENSION( n ), INTENT( IN ) :: X
      REAL ( c_double ), DIMENSION( m ), INTENT( OUT ) :: F
      class( params_base_type ), intent(in) :: params
      REAL ( c_double ) :: obj

!  evaluate the residuals F

      CALL CUTEST_cfn( status, n, m, X, obj, F )
      RETURN
      END SUBROUTINE eval_F

      SUBROUTINE eval_J( status, n, m, X, J, params)
      USE ISO_C_BINDING
      use :: nlls_module, only : params_base_type
      
      INTEGER ( c_int ), INTENT( OUT ) :: status
      INTEGER ( c_int ), INTENT( IN ) :: n, m
      REAL ( c_double ), DIMENSION( n ), INTENT( IN ) :: X
      REAL ( c_double ), DIMENSION( m * n ), INTENT( OUT ) :: J
      class( params_base_type ), intent(in) :: params
      REAL ( c_double ), DIMENSION( n ) :: G
      REAL ( c_double ), DIMENSION( m ) :: Y
      REAL ( c_double ), DIMENSION( m , n ) :: Jmatrix

!  evaluate the residual Jacobian J

      CALL CUTEST_cgr( status, n, m, X, Y, .FALSE., G, .FALSE., m, n, Jmatrix ) 
      ! convert the Jacobian to a vector....
      J = reshape(Jmatrix, (/n*m/) )
      RETURN
      END SUBROUTINE eval_J

      SUBROUTINE eval_HF( status, n, m, X, F, H, params)
      USE ISO_C_BINDING
      use :: nlls_module, only : params_base_type
      
      INTEGER ( c_int ), INTENT( OUT ) :: status
      INTEGER ( c_int ), INTENT( IN ) :: n, m
      REAL ( c_double ), DIMENSION( n ), INTENT( IN ) :: X
      REAL ( c_double ), DIMENSION( m ), INTENT( IN ) :: F
      REAL ( c_double ), DIMENSION( n*n ), INTENT( OUT ) :: H
      class( params_base_type ), intent(in) :: params
      
      real ( c_double ), dimension(n,n) :: Hmatrix
!  evaluate the product H = sum F_i Hessian F_i

      CALL CUTEST_cdhc( status, n, m, X, F, n, Hmatrix )
      H = reshape(Hmatrix, (/n*n/) )
      RETURN
      END SUBROUTINE eval_HF




