!     ( Last modified on 29 Jan 2013 at 14:15:00 )

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
  
     TYPE, PUBLIC :: NLLS_control_type
       INTEGER :: error = 6
       INTEGER :: out = 6
       INTEGER :: print_level = 0
     END TYPE NLLS_control_type

     TYPE, PUBLIC :: NLLS_inform_type
       INTEGER :: status = 0
       REAL( c_double ) :: obj = HUGE( 1.0_c_double )
     END TYPE NLLS_inform_type

   CONTAINS

     SUBROUTINE RAL_NLLS( n, m, X, Work_int, len_work_int,                     &
                          Work_real, len_work_real,                            &
                          eval_F, eval_J, eval_HF,                             &
                          control, inform )
    
!  -----------------------------------------------------------------------------
!  RAL_NLLS, a fortran subroutine for finding a first-order critical
!   point (most likely, a local minimizer) of the nonlinear least-squares 
!   objective function 1/2 ||F(x)||_2^2.

!  Authors: RAL NA Group (Iain Duff, Nick Gould, Jonathan Hogg, Tyrone Rees, 
!                         Jennifer Scott)
!  -----------------------------------------------------------------------------

!   Dummy arguments

     INTEGER( c_int ), INTENT( IN ) :: n, m, len_work_int, len_work_real
     REAL( c_double ), DIMENSION( n ), INTENT( INOUT ) :: X
     INTEGER( c_int ), INTENT( OUT ) :: Work_int( len_work_int )
     REAL( c_double ), INTENT( OUT ) :: Work_real( len_work_real )
     TYPE( NLLS_inform_type ), INTENT( OUT ) :: inform
     TYPE( NLLS_control_type ), INTENT( IN ) :: control

!  Interface blocks

     INTERFACE
       SUBROUTINE eval_F( status, n, m, X, F )
         USE ISO_C_BINDING
         INTEGER ( c_int ), INTENT( OUT ) :: status
         INTEGER ( c_int ), INTENT( IN ) :: n, m
         REAL ( c_double ), DIMENSION( n ), INTENT( IN ) :: X
         REAL ( c_double ), DIMENSION( m ), INTENT( OUT ) :: F
       END SUBROUTINE eval_F
     END INTERFACE

     INTERFACE
       SUBROUTINE eval_J( status, n, m, X, J )
         USE ISO_C_BINDING
         INTEGER ( c_int ), INTENT( OUT ) :: status
         INTEGER ( c_int ), INTENT( IN ) :: n, m
         REAL ( c_double ), DIMENSION( n ), INTENT( IN ) :: X
         REAL ( c_double ), DIMENSION( m , n ), INTENT( OUT ) :: J
       END SUBROUTINE eval_J
     END INTERFACE

     INTERFACE
       SUBROUTINE eval_HF( status, n, m, X, F, H )
         USE ISO_C_BINDING
         INTEGER ( c_int ), INTENT( OUT ) :: status
         INTEGER ( c_int ), INTENT( IN ) :: n, m
         REAL ( c_double ), DIMENSION( n ), INTENT( IN ) :: X
         REAL ( c_double ), DIMENSION( m ), INTENT( IN ) :: F
         REAL ( c_double ), DIMENSION( n , n ), INTENT( OUT ) :: H
       END SUBROUTINE eval_HF
     END INTERFACE

!  Local variables

     INTEGER :: status, start_f, end_f, start_j, start_h, w_end

!  check input dimensions

     IF ( m <= 0 .OR. n <= 0 ) THEN
       status = error_dimensions
       GO TO 990
     END IF

!  partition the workspace

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

     CALL eval_F( status, n, m, X, WORK_real( start_f ) )
     IF ( status /= 0 ) THEN
       status = error_eval_F
       GO TO 990
     END IF
     inform%obj = 0.5_c_double * DOT_PRODUCT( WORK_real( start_f : end_f ),    &
                                              WORK_real( start_f : end_f ) )

!  evaluate J

     CALL eval_J( status, n, m, X, WORK_real( start_j ) )
     IF ( status /= 0 ) THEN
       status = error_eval_J
       GO TO 990
     END IF

!  evaluate HF

     CALL eval_HF( status, n, m, X, WORK_real( start_f ), WORK_real( start_h ) )
     IF ( status /= 0 ) THEN
       status = error_eval_HF
       GO TO 990
     END IF

 990 CONTINUE
     inform%status = status
     RETURN
     END SUBROUTINE RAL_NLLS

   END MODULE NLLS_MODULE
