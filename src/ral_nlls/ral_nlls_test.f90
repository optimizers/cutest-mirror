!  Dummy RAL_NLLS for testing ral_nlls_main interface to CUTEst
!  Nick Gould, 6th October 2015

   MODULE NLLS_MODULE

     USE ISO_C_BINDING

     IMPLICIT none

     INTEGER, PARAMETER :: wp = KIND( 1.0d0 )
     INTEGER, PARAMETER :: long = SELECTED_INT_KIND( 8 )

     INTEGER, PARAMETER :: error_dimensions = - 1
     INTEGER, PARAMETER :: error_workspace = - 2
     INTEGER, PARAMETER :: error_eval_F = - 3
     INTEGER, PARAMETER :: error_eval_J = - 4
     INTEGER, PARAMETER :: error_eval_HF = - 5

     real (kind = wp), parameter :: tenm5 = 1.0e-5
     real (kind = wp), parameter :: tenm8 = 1.0e-8
     real (kind = wp), parameter :: epsmch = epsilon(1.0_wp)
     real (kind = wp), parameter :: hundred = 100.0
     real (kind = wp), parameter :: ten = 10.0
     real (kind = wp), parameter :: point9 = 0.9
     real (kind = wp), parameter :: zero = 0.0
     real (kind = wp), parameter :: one = 1.0
     real (kind = wp), parameter :: two = 2.0
     real (kind = wp), parameter :: half = 0.5
     real (kind = wp), parameter :: sixteenth = 0.0625
  
     TYPE, PUBLIC :: NLLS_control_type
       INTEGER :: error = 6
       INTEGER :: out = 6
       INTEGER :: print_level = 0
       INTEGER :: maxit = 100
       INTEGER :: model = 1
       INTEGER :: nlls_method = 1
       INTEGER :: lls_solver = 1
       REAL ( KIND = wp ) :: stop_g_absolute = tenm5
       REAL ( KIND = wp ) :: stop_g_relative = tenm8
       REAL ( KIND = wp ) :: initial_radius = hundred
       REAL ( KIND = wp ) :: maximum_radius = ten ** 8
       REAL ( KIND = wp ) :: eta_successful = ten ** ( - 8 )
       REAL ( KIND = wp ) :: eta_very_successful = point9
       REAL ( KIND = wp ) :: eta_too_successful = two
       REAL ( KIND = wp ) :: radius_increase = two
       REAL ( KIND = wp ) :: radius_reduce = half
       REAL ( KIND = wp ) :: radius_reduce_max = sixteenth

     LOGICAL :: subproblem_eig_fact = .FALSE.
     integer  :: more_sorensen_maxits = 500
     real(wp) :: more_sorensen_shift = 1e-8
     real(wp) :: more_sorensen_tiny = 10.0 * epsmch
     real(wp) :: more_sorensen_tol = 1e-6

     END TYPE NLLS_control_type

     TYPE, PUBLIC :: NLLS_inform_type
     INTEGER :: status = 0
     INTEGER :: alloc_status = 0
     INTEGER :: iter = 0
     INTEGER :: f_eval = 0
     INTEGER :: g_eval = 0
     INTEGER :: h_eval = 0
     REAL ( KIND = wp ) :: obj = HUGE( one )
     REAL ( KIND = wp ) :: norm_g = HUGE( one )
!      REAL( c_double ) :: obj = HUGE( 1.0_c_double )
     END TYPE NLLS_inform_type

     type params_base_type
     ! deliberately empty
     end type params_base_type

   CONTAINS

     SUBROUTINE RAL_NLLS( n, m, X, eval_F, eval_J, eval_HF,                    &
                          params, inform, control )
    
!  -----------------------------------------------------------------------------
!  RAL_NLLS, a fortran subroutine for finding a first-order critical
!   point (most likely, a local minimizer) of the nonlinear least-squares 
!   objective function 1/2 ||F(x)||_2^2.

!  Authors: RAL NA Group (Iain Duff, Nick Gould, Jonathan Hogg, Tyrone Rees, 
!                         Jennifer Scott)
!  -----------------------------------------------------------------------------

!   Dummy arguments

     INTEGER( c_int ), INTENT( IN ) :: n, m
     REAL( c_double ), DIMENSION( n ), INTENT( INOUT ) :: X
     class( params_base_type ) :: params
     TYPE( NLLS_inform_type ), INTENT( OUT ) :: inform
     TYPE( NLLS_control_type ), INTENT( IN ) :: control

!  Interface blocks

     INTERFACE
       SUBROUTINE eval_F( status, n, m, X, F , params )
         USE ISO_C_BINDING
         import :: params_base_type
         INTEGER ( c_int ), INTENT( OUT ) :: status
         INTEGER ( c_int ), INTENT( IN ) :: n, m
         REAL ( c_double ), DIMENSION( n ), INTENT( IN ) :: X
         REAL ( c_double ), DIMENSION( m ), INTENT( OUT ) :: F
         class( params_base_type ), intent( in ) :: params
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
         class( params_base_type ), intent( in ) :: params
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
         REAL ( c_double ), DIMENSION( n*n ), INTENT( OUT ) :: H
         class( params_base_type ), intent( in ) :: params
       END SUBROUTINE eval_HF
     END INTERFACE

!  Local variables

     INTEGER :: status, start_f, end_f, start_j, start_h, w_end
     INTEGER :: len_work_int, len_work_real
     INTEGER( c_int ), allocatable :: Work_int( : )
     REAL( c_double ), allocatable :: Work_real( : ) 
     
!  check input dimensions

     IF ( m <= 0 .OR. n <= 0 ) THEN
       status = error_dimensions
       GO TO 990
     END IF

!  partition the workspace
     allocate(Work_int(10))
     allocate(Work_real(m + n*(n + m)))
     
     start_f = 1
     start_j = start_f + m
     end_f = start_j - 1
     start_h = start_j + n * m
     w_end = start_h + n * n - 1

     IF ( w_end < len_work_real ) THEN
       status = error_workspace
       GO TO 990
     END IF

!  evaluate F

     CALL eval_F( status, n, m, X, WORK_real( start_f ), params )
     IF ( status /= 0 ) THEN
       status = error_eval_F
       GO TO 990
     END IF
     inform%obj = 0.5_c_double * DOT_PRODUCT( WORK_real( start_f : end_f ),    &
                                              WORK_real( start_f : end_f ) )

!  evaluate J

     CALL eval_J( status, n, m, X, WORK_real( start_j ), params )
     IF ( status /= 0 ) THEN
       status = error_eval_J
       GO TO 990
     END IF

!  evaluate HF

     CALL eval_HF( status, n, m, X, WORK_real( start_f ), WORK_real( start_h ), params )
     IF ( status /= 0 ) THEN
       status = error_eval_HF
       GO TO 990
     END IF

 990 CONTINUE
     inform%status = status
     RETURN
     END SUBROUTINE RAL_NLLS

   END MODULE NLLS_MODULE
