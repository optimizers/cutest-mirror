! THIS VERSION: CUTEST 1.1 - 16/01/2014 AT 10:15 GMT.

!-*-*-*-*-*-*- C U T E S T   q p l i b  _ m a i n   P R O G R A M -*-*-*-*-*-

    PROGRAM CUTEST_qplib_main

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 2003 version released January 2014

     USE CUTEST_LQP_double

!--------------------
!   P r e c i s i o n
!--------------------

      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!----------------------
!   P a r a m e t e r s
!----------------------

      INTEGER, PARAMETER :: input = 55
      INTEGER, PARAMETER :: input_spec = 46
      INTEGER, PARAMETER :: standard_out = 6
      INTEGER, PARAMETER :: qplib_out = 61
      INTEGER, PARAMETER :: buffer = 77 
      REAL ( KIND = wp ), PARAMETER :: zero = 0.0_wp
      REAL ( KIND = wp ), PARAMETER :: one = 1.0_wp
      REAL ( KIND = wp ), PARAMETER :: infinity = ( 10.0_wp ) ** 19

      INTEGER, PARAMETER :: qp = 1
      INTEGER, PARAMETER :: qpqc = 2
      INTEGER, PARAMETER :: bqp = 3
      INTEGER, PARAMETER :: lp = 4
      INTEGER, PARAMETER :: lpqc = 5

!     CHARACTER ( len = 16 ) :: char_int_default = REPEAT( ' ', 16 )
!     CHARACTER ( len = 24 ) :: char_val_default = REPEAT( ' ', 24 )

!--------------------------------
!   L o c a l   V a r i a b l e s
!--------------------------------

      INTEGER :: n, m, H_ne, A_ne, status, out
      REAL ( KIND = wp ) :: f, h_pert
      CHARACTER ( len = 10 ) ::  p_name
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: X_type
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: A_row, A_col
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: H_row, H_col
      REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: X, X_l, X_u, Z, G
      REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: Y, C_l, C_u
      REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: A_val, H_val
      CHARACTER ( len = 10 ), ALLOCATABLE, DIMENSION( : ) :: X_names, C_names

      INTEGER :: i, int_var, l
      INTEGER :: problem_type = 1
      LOGICAL :: filexx
      LOGICAL :: append_dim = .FALSE.
      LOGICAL :: qplib_wrfile = .FALSE.
      CHARACTER ( len = 16 ) :: char_i, char_j, char_l
      CHARACTER ( len = 24 ) :: char_val
      CHARACTER ( len = 28 ) :: out_p_name
      CHARACTER ( len = 34 ) :: out_p_name_qplib
      REAL ( KIND = wp ) :: mode_v
!     REAL ( KIND = wp ), DIMENSION( 100 ) :: V

      out = standard_out

!  open the Spec file for the method

      OPEN( input_spec, FILE = 'QPLIB.SPC', FORM = 'FORMATTED', STATUS = 'OLD')
      REWIND( input_spec )

!  read input Spec data

!  problem_type = 1(QP,default),2,(QPQC),3(BQP),4(LP),5(LPQC)
!  append_dim   = F(don't append _n_m to name where n.m = dims,default), T(do)
!  qplib_wrfile = T(write to a file "probname".qplib), F(write to standard out) 

!  set up algorithmic input data

      READ ( input_spec, * ) problem_type
      READ ( input_spec, * ) append_dim
      READ ( input_spec, * ) qplib_wrfile
      CLOSE ( input_spec )

!  open the problem data file

      OPEN ( input, FILE = 'OUTSDIF.d', FORM = 'FORMATTED', STATUS = 'OLD' )

!  sparse co-ordinate version

      h_pert = zero
      IF ( problem_type == qpqc .OR. problem_type == lpqc ) THEN
        CALL CUTEST_lqp_create( status, input, buffer, out, n, m, f, G,        &
                                X, X_l, X_u, Z, Y, C_l, C_u, p_name,           &
                                X_names, C_names,                              &
                                X_type = X_type,                               &
                                A_ne = A_ne, A_row = A_row,                    &
                                A_col = A_col, A_val = A_val,                  &
                                H_ne = H_ne, H_row = H_row,                    &
                                H_col = H_col, H_val = H_val,                  &
                                H_pert = h_pert )
      ELSE
        CALL CUTEST_lqp_create( status, input, buffer, out, n, m, f, G,        &
                                X, X_l, X_u, Z, Y, C_l, C_u, p_name,           &
                                X_names, C_names,                              &
                                X_type = X_type,                               &
                                A_ne = A_ne, A_row = A_row,                    &
                                A_col = A_col, A_val = A_val,                  &
                                H_ne = H_ne, H_row = H_row,                    &
                                H_col = H_col, H_val = H_val,                  &
                                H_pert = h_pert )
      END IF
      CLOSE( input )
      IF ( status /= 0 ) GO TO 900
!      CALL WRITE_p_name( out, p_name )
!      WRITE( out, "( ' * n = ', I0, ', m = ', I0,                             &
!     &  ', A_ne = ', I0, ', H_ne = ', I0 )" ) n, m, A_ne, H_ne
!     CALL WRITE_H_sparse( out, H_ne, H_ne, H_val, H_row, H_col )
!     CALL WRITE_G( out, n, G )
!     CALL WRITE_f( out, f )
!     CALL WRITE_A_sparse( out, A_ne, A_ne, A_val, A_row, A_col )
!     CALL WRITE_X( out, n, X, X_l, X_u, Z )
!     CALL WRITE_Y( out, m, Y, C_l, C_u )
!     CALL WRITE_X_type( out, n, X_type )
!     CALL WRITE_X_names( out, n, X_names )
!     CALL WRITE_C_names( out, m, C_names )

!  obtain output problem name

      out_p_name = REPEAT( ' ', 28 )
      IF ( append_dim ) THEN
        char_l = TRIM_INT( n )
        IF ( m > 0 ) THEN
          char_i = TRIM_INT( m )
          out_p_name = TRIM( p_name ) // '_' // TRIM( char_l )                 &
                                      // '_' // TRIM( char_i )
        ELSE
          out_p_name = TRIM( p_name ) // '_' // TRIM( char_l )
        END IF
      ELSE
        out_p_name = TRIM( p_name )
      END IF

!  open output file if required

      IF ( qplib_wrfile ) THEN
        out_p_name_qplib = REPEAT( ' ', 34 )
        out_p_name_qplib = TRIM( out_p_name ) // ".qplib"
        INQUIRE( FILE = out_p_name_qplib, EXIST = filexx )
        IF ( filexx ) THEN
           OPEN( qplib_out, FILE = out_p_name_qplib, FORM = 'FORMATTED',       &
                 STATUS = 'OLD', IOSTAT = status )
        ELSE
           OPEN( qplib_out, FILE = out_p_name_qplib, FORM = 'FORMATTED',       &
                 STATUS = 'NEW', IOSTAT = status )
        END IF
        IF ( status /= 0 ) GO TO 900
        out = qplib_out
        REWIND ( out )
      END IF

!  see if the problem has integer variables

      int_var = COUNT( X_type( : n ) > 0 )

!  set header

      IF ( m > 0 ) THEN
        IF ( LEN( TRIM( out_p_name ) ) <= 24 ) THEN
          WRITE( out, "( A24, ' CUTEst ', A, '.SIF with n = ', I0, ', m = ',   &
         &   I0 )" ) out_p_name( 1 : 24 ), TRIM( p_name ), n, m
        ELSE
          WRITE( out, "( A28, /, 24X, ' CUTEst ', A, '.SIF with n = ', I0,     &
         & ', m = ', I0 )" ) out_p_name, TRIM( p_name ), n, m
        END IF
      ELSE
        WRITE( out, "( A24, ' CUTEst ', A,                                     &
       &       '.SIF with n = ', I0 )" ) out_p_name, TRIM( p_name ), n
      END IF 
      IF ( int_var == 0 ) THEN
        SELECT CASE ( problem_type )
        CASE ( qpqc )
          WRITE( out, "( 'QPQC                     a quadratic program',       &
         &               ' with quadratic constraints' )" )
        CASE ( bqp )
          WRITE( out, "( 'BQP                      a bound-constrained',       &
         &               ' quadratic program' )" )
        CASE ( lp )
          WRITE( out, "( 'LP                       a linear program' )" )
        CASE ( lpqc )
          WRITE( out, "( 'LPQC                     a linear program',          &
         &                ' with quadratic constraints' )" )
        CASE DEFAULT
          WRITE( out, "( 'QP                       a quadratic program' )")
        END SELECT
      ELSE IF ( int_var == n ) THEN
        SELECT CASE ( problem_type )
        CASE ( qpqc )
          WRITE( out, "( 'IQPQC                    an integer',                &
         &    ' QP with quadratic constraints' )" )
        CASE ( bqp )
          WRITE( out, "( 'IBQP                     an integer',                &
         &    ' bound-constrained quadratic program' )" )
        CASE ( lp )
          WRITE( out, "( 'ILP                      an integer',                &
       &     ' linear program' )" )
        CASE ( lpqc )
          WRITE( out, "( 'ILPQC                    an integer',                &
         &    ' LP with quadratic constraints' )" )
        CASE DEFAULT
          WRITE( out, "( 'IQP                      an integer',                &
         &   ' quadratic program' )")
        END SELECT
      ELSE
        SELECT CASE ( problem_type )
        CASE ( qpqc )
          WRITE( out, "( 'MIQPQC                   a mixed-integer',           &
         &    ' QP with quadratic constraints' )" )
        CASE ( bqp )
          WRITE( out, "( 'MIBQP                    a mixed-integer',           &
         &    ' bound-constrained quadratic program' )" )
        CASE ( lp )
          WRITE( out, "( 'MILP                     a mixed-integer',           &
       &     ' linear program' )" )
        CASE ( lpqc )
          WRITE( out, "( 'MILPQC                   a mixed-integer',           &
         &    ' LP with quadratic constraints' )" )
        CASE DEFAULT
          WRITE( out, "( 'MIQP                     a mixed-integer',           &
         &   ' quadratic program' )")
        END SELECT
      END IF
      char_l = TRIM_INT( n )
      WRITE( out, "( A16, 8X, ' # variables ' )" ) char_l
      IF ( problem_type /= bqp ) THEN
        char_l = TRIM_INT( m )
        WRITE( out, "( A16, 8X, ' # general linear constraints ' )" ) char_l
      END IF

!  Hessian values

      IF ( problem_type == qp .OR. problem_type == bqp .OR.                    &
           problem_type == qpqc ) THEN
        char_l = TRIM_INT( H_ne )
        IF ( H_ne == 0 ) THEN
          WRITE( out, "( /, A16, 8X, ' # nonzeros in upper triangle of H' )" ) &
            char_l
        ELSE
          WRITE( out, "( /, A16, 8X, ' # nonzeros in upper triangle of H:',    &
         &   ' row,column,value' )" ) char_l
          DO l = 1, H_ne
            char_i = TRIM_INT( H_row( l ) ) ; char_j = TRIM_INT( H_col( l ) )
            char_val = TRIM_VALUE( H_val( l ) )
            WRITE( out, "( A16, 1X, A16, 1X, A24 ) )" ) char_i, char_j, char_val
          END DO
        END IF
      END IF

!  gradient values

      mode_v = MODE( n, G )
      l = COUNT( G( : n ) /= mode_v )
      char_l = TRIM_INT( l ) ; char_val = TRIM_VALUE( mode_v )
      WRITE( out, "( /, A24, ' default value for entries in g' )" ) char_val
      IF ( l == 0 ) THEN
        WRITE( out, "( A16, 8X, ' # non default entries in g' )" ) char_l
      ELSE
        WRITE( out, "( A16, 8X, ' # non default entries in g: index,value' )") &
          char_l
        DO i = 1, n
          IF ( G( i ) /= mode_v ) THEN
            char_i = TRIM_INT( i ) ; char_val = TRIM_VALUE( G( i ) )
            WRITE( out, "( A16, 1X, A24 ) )" ) char_i, char_val
          END IF
        END DO
      END IF

!  function value

      char_val = TRIM_VALUE( f )
      WRITE( out, "( /, A24, ' value of f' ) )" ) char_val

!  Hessian values for constraints

      IF ( problem_type == qpqc .OR. problem_type == lpqc ) THEN
        char_l = TRIM_INT( 0 )
        WRITE( out, "( /, A16, 8X, ' # nonzeros in upper triangle',            &
       &  ' of the H_i')") char_l
      END IF

!  constraint Jacobian values

      IF ( problem_type /= bqp ) THEN
        char_l = TRIM_INT( A_ne )
        IF ( A_ne == 0 ) THEN
          WRITE( out, "( /, A16, 8X, ' # nonzeros in A' )" ) char_l
        ELSE
          WRITE( out, "( /, A16, 8X, ' # nonzeros in A:',                      &
         &   ' row,column,value' )" ) char_l
          DO l = 1, A_ne
            char_i = TRIM_INT( A_row( l ) ) ; char_j = TRIM_INT( A_col( l ) )
            char_val = TRIM_VALUE( A_val( l ) )
            WRITE( out, "( A16, 1X, A16, 1X, A24 ) )" ) char_i, char_j, char_val
          END DO
        END IF
      END IF

!  infinity

      char_val = TRIM_VALUE( infinity )
      WRITE( out, "( /, A24, ' value of infinite bounds' ) )" ) char_val

!  constraint lower bounds
    
      IF ( problem_type /= bqp ) THEN
        mode_v = MODE( m, C_l )
        l = COUNT( C_l( : m ) /= mode_v )
        char_l = TRIM_INT( l ) ; char_val = TRIM_VALUE( mode_v )
        WRITE( out, "( /, A24, ' default value for entries in c_l' )" ) char_val
        IF ( l == 0 ) THEN
          WRITE( out, "( A16, 8X, ' # non default entries in c_l' )" ) char_l
        ELSE
          WRITE( out, "( A16, 8X, ' # non default entries in c_l:',            &
         &  ' index,value' )" ) char_l
          DO i = 1, m
            IF ( C_l( i ) /= mode_v ) THEN
              char_i = TRIM_INT( i ) ; char_val = TRIM_VALUE( C_l( i ) )
              WRITE( out, "( A16, 1X, A24 ) )" ) char_i, char_val
            END IF
          END DO
        END IF

!  constraint upper bounds
    
        mode_v = MODE( m, C_u )
        l = COUNT( C_u( : m ) /= mode_v )
        char_l = TRIM_INT( l ) ; char_val = TRIM_VALUE( mode_v )
        WRITE( out, "( /, A24, ' default value for entries in c_u' )" ) char_val
        IF ( l == 0 ) THEN
          WRITE( out, "( A16, 8X, ' # non default entries in c_u' )" ) char_l
        ELSE
          WRITE( out, "( A16, 8X, ' # non default entries in c_u:',            &
        &   ' index,value' )" ) char_l
          DO i = 1, m
            IF ( C_u( i ) /= mode_v ) THEN
              char_i = TRIM_INT( i ) ; char_val = TRIM_VALUE( C_u( i ) )
              WRITE( out, "( A16, 1X, A24 ) )" ) char_i, char_val
            END IF
          END DO
        END IF
      END IF

!  variable lower bounds
    
      mode_v = MODE( n, X_l )
      l = COUNT( X_l( : n ) /= mode_v )
      char_l = TRIM_INT( l ) ; char_val = TRIM_VALUE( mode_v )
      WRITE( out, "( /, A24, ' default value for entries in x_l' )" ) char_val
      IF ( l == 0 ) THEN
        WRITE( out, "( A16, 8X, ' # non default entries in x_l' )" ) char_l
      ELSE
        WRITE( out, "( A16, 8X, ' # non default entries in x_l:',              &
       &  ' index,value' )" ) char_l
        DO i = 1, n
          IF ( X_l( i ) /= mode_v ) THEN
            char_i = TRIM_INT( i ) ; char_val = TRIM_VALUE( X_l( i ) )
            WRITE( out, "( A16, 1X, A24 ) )" ) char_i, char_val
          END IF
        END DO
      END IF

!  variable upper bounds
    
      mode_v = MODE( n, X_u )
      l = COUNT( X_u( : n ) /= mode_v )
      char_l = TRIM_INT( l ) ; char_val = TRIM_VALUE( mode_v )
      WRITE( out, "( /, A24, ' default value for entries in x_u' )" ) char_val
      IF ( l == 0 ) THEN
        WRITE( out, "( A16, 8X, ' # non default entries in x_u' )" ) char_l
      ELSE
        WRITE( out, "( A16, 8X, ' # non default entries in x_u:',              &
       &  ' index,value' )") char_l
        DO i = 1, n
          IF ( X_u( i ) /= mode_v ) THEN
            char_i = TRIM_INT( i ) ; char_val = TRIM_VALUE( X_u( i ) )
            WRITE( out, "( A16, 1X, A24 ) )" ) char_i, char_val
          END IF
        END DO
      END IF

!  variable types

      IF ( int_var > 0 .AND. int_var < n ) THEN
        IF ( n >= 2 * int_var ) THEN
          char_l = TRIM_INT( 0 )
          WRITE( out, "( /, A16, 8X, ' default variable type',                 &
         &  ' (0 for continuous, 1 for integer)' )" ) char_l
          char_j = TRIM_INT( int_var )
          IF ( int_var == 0 ) THEN
            WRITE( out, "( A16, 8X, ' # non default variables' )" ) char_j
          ELSE
            WRITE( out, "( A16, 8X, ' # non default variables: index,type' )") &
           &  char_j
            DO i = 1, n
              IF (  X_type( i ) /= 0 ) THEN
                char_i = TRIM_INT( i ) ; char_j = TRIM_INT( X_type( i ) )
                WRITE( out, "( A16, 1X, A16 ) )" ) char_i, char_j
              END IF
            END DO
          END IF
        ELSE
          char_l = TRIM_INT( 1 )
          WRITE( out, "( /, A16, 8X, ' default variable type',                 &
         & ' (0 for continuous, 1 for integer)' )" ) char_l
          char_j = TRIM_INT( n - int_var )
          IF ( int_var == n ) THEN
            WRITE( out, "( A16, 8X, ' # non default variables' )" ) char_j
          ELSE
            WRITE( out, "( A16, 8X, ' # non default variables: index,type' )") &
              char_j
            DO i = 1, n
              IF (  X_type( i ) == 0 ) THEN
                char_i = TRIM_INT( i ) ; char_j = TRIM_INT( X_type( i ) )
                WRITE( out, "( A16, 1X, A16 ) )" ) char_i, char_j
              END IF
            END DO
          END IF
        END IF
      END IF

!  initial primal variables

      mode_v = MODE( n, X )
      l = COUNT( X( : n ) /= mode_v )
      char_l = TRIM_INT( l ) ; char_val = TRIM_VALUE( mode_v )
      WRITE( out, "( /, A24, ' default value for entries in initial x' )" )    &
        char_val
      IF ( l == 0 ) THEN
        WRITE( out, "( A16, 8X, ' # non default entries in x' )" ) char_l
      ELSE
        WRITE( out, "( A16, 8X, ' # non default entries in x: index,value' )") &
          char_l
        DO i = 1, n
          IF ( X( i ) /= mode_v ) THEN
            char_i = TRIM_INT( i ) ; char_val = TRIM_VALUE( X( i ) )
            WRITE( out, "( A16, 1X, A24 ) )" ) char_i, char_val
          END IF
        END DO
      END IF

!  initial Lagrange multipliers

      IF ( problem_type /= bqp ) THEN
        mode_v = MODE( m, Y )
        l = COUNT( Y( : m ) /= mode_v )
        char_l = TRIM_INT( l ) ; char_val = TRIM_VALUE( mode_v )
        WRITE( out, "( /, A24, ' default value for entries in initial y' )" )  &
          char_val
        IF ( l == 0 ) THEN
          WRITE( out, "( A16, 8X, ' # non default entries in y' )" ) char_l
        ELSE
          WRITE( out, "( A16, 8X, ' # non default entries in y:',              &
         &  ' index,value' )" ) char_l
          DO i = 1, m
            IF ( Y( i ) /= mode_v ) THEN
              char_i = TRIM_INT( i ) ; char_val = TRIM_VALUE( Y( i ) )
              WRITE( out, "( A16, 1X, A24 ) )" ) char_i, char_val
            END IF
          END DO
        END IF
      END IF

!  initial dual variables

      mode_v = MODE( n, Z )
      l = COUNT( Z( : n ) /= mode_v )
      char_l = TRIM_INT( l ) ; char_val = TRIM_VALUE( mode_v )
      WRITE( out, "( /, A24, ' default value for entries in initial z' )" )    &
        char_val
      IF ( l == 0 ) THEN
        WRITE( out, "( A16, 8X, ' # non default entries in z' )" ) char_l
      ELSE
        WRITE( out, "( A16, 8X, ' # non default entries in z: index,value' )") &
          char_l
        DO i = 1, n
          IF ( Z( i ) /= mode_v ) THEN
            char_i = TRIM_INT( i ) ; char_val = TRIM_VALUE( Z( i ) )
            WRITE( out, "( A16, 1X, A24 ) )" ) char_i, char_val
          END IF
        END DO
      END IF

!  variable names

      char_l = TRIM_INT( n )
      WRITE( out, "( /, A16, 8X, ' # non default names for variables:',        &
     &   ' index,name' )" ) char_l
      DO i = 1, n
        char_i = TRIM_INT( i )
        WRITE( out, "( A16, 1X, A10 ) )" ) char_i, X_names( i )
      END DO

!  constraint names

      IF ( problem_type /= bqp ) THEN
        char_l = TRIM_INT( n )
        WRITE( out, "( /, A16, 8X, ' # non default names for constraints:',    &
       &   ' index,name' )" ) char_l
        DO i = 1, m
          char_i = TRIM_INT( i )
          WRITE( out, "( A16, 1X, A10 ) )" ) char_i, C_names( i )
        END DO
      END IF

      DEALLOCATE( A_row, A_col, H_row, H_col, X, X_l, X_u, Z, G, Y, C_l, C_u,  &
                  A_val, H_val, X_names, C_names, X_type, stat = status )

      IF ( qplib_wrfile ) CLOSE( qplib_out )

      STOP

!  error exits

 900  CONTINUE
      WRITE( out, "( ' error status = ', I0 )" ) status
      CLOSE( INPUT  )
      STOP

    CONTAINS

      FUNCTION MODE( n, V )
      IMPLICIT NONE
      REAL ( KIND = wp ) :: MODE
      INTEGER, INTENT( IN ) :: n
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: V

!  find the "mode", i.e., the most commonly-occuring value, of a vector v

      INTEGER :: i, mode_start, max_len, same, len, m, inform

      REAL ( KIND = wp ), DIMENSION( n ) :: V_sorted

!  sort a copy of v into increasing order

      V_sorted = V
      CALL SORT_heapsort_build( n, V_sorted, inform ) !  build the heap
      DO i = 1, n
        m = n - i + 1 
        CALL SORT_heapsort_smallest( m, V_sorted, inform ) !  reorder v
      END DO  

!  run through the sorted values, finding adjacent entries that are identical

      mode_start = 1 ; max_len = 1
      same = 1 ; len = 1
      DO i = 2, n
        IF ( V_sorted( i ) /= V_sorted( same ) ) THEN
          IF ( len > max_len ) THEN
            mode_start = same
            max_len = len
          END IF
          same = i ; len = 1
        ELSE
          len = len + 1
        END IF
      END DO
      IF ( len > max_len ) THEN
        mode_start = same
        max_len = len
      END IF
!     write(6,*) max_len
!     write(6,*) V_sorted( : n )
      MODE = V_sorted( mode_start )
      RETURN

      END FUNCTION MODE

!  heapsort routines extracted from GALAHAD

      SUBROUTINE SORT_heapsort_build( n, A, inform )

!  Given an array A, elements A(1), ...., A(N), subroutine SORT_heapsort_build
!  rearranges the elements to form a heap in which each parent has a smaller
!  value than either of its children.

!  Algorithm 232 of CACM (J. W. J. Williams): a combination of
!  ALGOL procedures SETHEAP and INHEAP

!  Programming: Nick Gould, January 26th 1995.

!  ------------------------- dummy arguments --------------------------
!
!  n      integer, which gives the number of values to be sorted.
!         n must be positive
!
!  A      real array of length n. On input, A must contain the
!         values which are to be sorted. On output, these values
!         will have been permuted so as to form a heap
!
!  inform integer, which informs the user of the success of SORT_heapsort_build.
!         If inform = 0 on exit, the heap has been formed.
!         If inform = 1 on exit, n was input with a value less than
!                       or equal to 0 and the heap has not been formed.
!
!  ------------------ end of dummy arguments --------------------------

      IMPLICIT NONE
      INTEGER, INTENT( IN ) :: n
      INTEGER, INTENT( OUT ) :: inform
      REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( n ) :: A

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: i, j, k
      REAL ( KIND = wp )  :: rin

!  Add the elements to the heap one at a time

      IF ( n <= 0 ) THEN
         inform = 1
         RETURN
      ENDIF

      DO k = 2, n
        rin = A( k )

!  The cycle may be repeated log2(k) times, but on average is repeated
!  merely twice

        i = k
        DO
          IF ( i <= 1 ) EXIT
          j = i / 2
          IF ( A( j ) <= rin ) EXIT
          A( i ) = A( j )
          i = j
        END DO
        A( i ) = rin
      END DO
      inform = 0

      RETURN

!  End of subroutine SORT_heapsort_build

     END SUBROUTINE SORT_heapsort_build

     SUBROUTINE SORT_heapsort_smallest( m, A, inform )

!  Given an array A, elements A(1), ...., A(m) forming a heap,
!  SORT_heapsort_smallest assigns to rout the value of A(1), the smallest
!  member of the heap, and arranges the remaining members as elements
!  1 to m - 1 of A. rout is then placed in A(m)

!  Algorithm 232 of CACM (J. W. J. Williams): a combination of
!  ALGOL procedures OUTHEAP and SWOPHEAP

!  Programming: Nick Gould, January 26th 1995.

!  ------------------------- dummy arguments --------------------------
!
!  m      integer, which gives the number of values to be sorted.
!         m must be positive
!
!  A      real array of length m. On input, A must contain the values which
!         are to be sorted stored in a heap. On output, the smallest value
!         will have been moved into A(m) and the remaining values A(k),
!         k = 1,..., m-1 will have been restored to a heap
!
!  inform integer, which informs the user of the success of 
!         SORT_heapsort_smallest.
!         If inform = 0 on exit, the smallest value has been found.
!         If inform = 1 on exit, m was input with a value less than
!                       or equal to 0 and the heap has not been formed
!
!  ------------------ end of dummy arguments --------------------------

      IMPLICIT NONE
      INTEGER, INTENT( IN ) :: m
      INTEGER, INTENT( OUT ) :: inform
      REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( m ) :: A

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: i, j
      REAL ( KIND = wp ) :: rin, rout

!  Add the element rin to the heap, extract and assign to rout
!  the value of the smallest member of the resulting set, and
!  leave the remaining elements in a heap of the original size.
!  In this process, elements 1 to n+1 of the array A may be disturbed

      IF ( m <= 0 ) THEN
         inform = 1
         RETURN
      ENDIF

      IF ( m > 1 ) THEN
        i = 1
        rout = A( 1 )
        rin = A( m )

!  Move from the top of the heap comparing the value of node i
!  with its two daughters. If node i is smallest, the heap has been
!  restored. If one of the children is smallest, promote this child
!  in the heap and move to the now vacated node.
!  This cycle may be repeated log2(m) times

        DO
          j = i + i
          IF ( j > m - 1 ) EXIT

!  Determine which of the two daughters is smallest

          IF ( A( j + 1 ) < A( j ) ) j = j + 1

!  Determine if the smaller daughter is less than the value from node i

          IF ( A( j ) >= rin ) EXIT
          A( i ) = A( j )
          i = j
        END DO

!  The heap has been restored

        A( i ) = rin

!  Store the smallest value in the now vacated m-th position of the list

        A( m ) = rout

      END IF
      inform = 0

      RETURN

!  End of subroutine SORT_heapsort_smallest

      END SUBROUTINE SORT_heapsort_smallest

      FUNCTION TRIM_VALUE( value )
      CHARACTER ( LEN = 24 ) :: TRIM_VALUE
      REAL ( KIND = wp ) :: value

!  remove more than one trailing 0 from floating point strings

      INTEGER :: i, j, k, zs
      LOGICAL :: zeros

      TRIM_VALUE = REPEAT( ' ', 24 )
      IF ( value > - 10.0_wp .AND. value < 10.0_wp ) THEN
        WRITE( TRIM_VALUE, "( F19.16 )" ) value
      ELSE
        WRITE( TRIM_VALUE, "( ES23.16 )" ) value
      END IF

!  remove any leading space

      IF ( TRIM_VALUE( 1 : 1 ) == ' ' ) THEN
        DO i = 2, 24
          TRIM_VALUE( i - 1 : i - 1 ) = TRIM_VALUE( i : i )
        END DO
      END IF 

      zeros = .FALSE.
      DO i = 1, 24
        IF ( TRIM_VALUE( i : i ) == '0' ) THEN
          IF ( .NOT. zeros ) THEN
            zs = i
            zeros = .TRUE.
          END IF
        ELSE IF ( TRIM_VALUE( i : i ) == 'E' .OR.                              &
                  TRIM_VALUE( i : i ) == 'e' .OR.                              &
                  TRIM_VALUE( i : i ) == 'D' .OR.                              &
                  TRIM_VALUE( i : i ) == 'd' ) THEN
          IF ( zeros ) THEN
            DO j = zs + 1, zs + 25 - i
              k = i + ( j - zs - 1 )
              TRIM_VALUE( j : j ) = TRIM_VALUE( k : k  )
            END DO
            DO j = zs + 26 - i, 24
              TRIM_VALUE( j : j ) = ' '
            END DO
          END IF
          zeros = .FALSE.
          EXIT
        ELSE IF ( TRIM_VALUE( i : i ) == ' ' ) THEN
          IF ( zeros ) THEN
            DO j = zs + 1, i
              TRIM_VALUE( j : j ) = ' '
            END DO
          END IF
          zeros = .FALSE.
          EXIT
        ELSE
          zeros = .FALSE.
        END IF
      END DO
      IF ( zeros ) THEN
        DO j = zs + 1, i
          TRIM_VALUE( j : j ) = ' '
        END DO
      END IF

!  remove superflous 0 from the exponent

      DO i = 1, 24
        IF ( TRIM_VALUE( i : i ) == 'E' .OR.                                   &
             TRIM_VALUE( i : i ) == 'e' .OR.                                   &
             TRIM_VALUE( i : i ) == 'D' .OR.                                   &
             TRIM_VALUE( i : i ) == 'd' ) THEN
          IF ( TRIM_VALUE( i + 1 : i + 1 ) == '+' .OR.                         &
               TRIM_VALUE( i + 1 : i + 1 ) == '-' ) THEN
            IF ( TRIM_VALUE( i + 2 : i + 2 ) == '0' ) THEN
              IF ( TRIM_VALUE( i + 3 : i + 3 ) == '0' ) THEN
                IF ( TRIM_VALUE( i + 4 : i + 4 ) == ' ' ) THEN
                  TRIM_VALUE( i + 3 : i + 3 ) = '0'
                ELSE
                  TRIM_VALUE( i + 2 : i + 2 ) = TRIM_VALUE( i + 4 : i + 4 )
                  TRIM_VALUE( i + 3 : i + 4 ) = '  '
                END IF
              ELSE
                IF ( TRIM_VALUE( i + 4 : i + 4 ) == ' ' ) THEN
                ELSE
                  TRIM_VALUE( i + 2 : i + 2 ) = TRIM_VALUE( i + 3 : i + 3 )
                  TRIM_VALUE( i + 3 : i + 3 ) = TRIM_VALUE( i + 4 : i + 4 )
                  TRIM_VALUE( i + 4 : i + 4 ) = ' '
                END IF
              END IF
            END IF
          ELSE
            IF ( TRIM_VALUE( i + 1 : i + 1 ) == '0' ) THEN
              IF ( TRIM_VALUE( i + 2 : i + 2 ) == '0' ) THEN
                IF ( TRIM_VALUE( i + 3 : i + 3 ) == ' ' ) THEN
                  TRIM_VALUE( i + 2 : i + 2 ) = '0'
                ELSE
                  TRIM_VALUE( i + 1 : i + 1 ) = TRIM_VALUE( i + 3 : i + 3 )
                  TRIM_VALUE( i + 2 : i + 3 ) = '  '
                END IF
              ELSE
                IF ( TRIM_VALUE( i + 3 : i + 3 ) == ' ' ) THEN
                ELSE
                  TRIM_VALUE( i + 1 : i + 1 ) = TRIM_VALUE( i + 2 : i + 2 )
                  TRIM_VALUE( i + 2 : i + 2 ) = TRIM_VALUE( i + 3 : i + 3 )
                  TRIM_VALUE( i + 3 : i + 3 ) = ' '
                END IF
              END IF
            END IF
          END IF
          EXIT
        END IF

!  remove trailing 0 unless it is preceeded by a .

        IF ( TRIM_VALUE( i : i ) == ' ' ) THEN
          IF ( i < 3 ) EXIT
          IF ( TRIM_VALUE( i - 1 : i - 1 ) == '0' .AND.                        &
               TRIM_VALUE( i - 2 : i - 2 ) /= '.' ) THEN
               TRIM_VALUE( i - 1 : i - 1 ) = ' '
          END IF
          EXIT
        END IF

      END DO

!  if the string starts with a ., add a 0 at the front

      IF ( TRIM_VALUE( 1 : 1 ) == '.' ) THEN
        DO i = 24, 2, -1
          TRIM_VALUE( i : i ) = TRIM_VALUE( i - 1 : i - 1 )
        END DO
        TRIM_VALUE( 1 : 1 ) = '0'
      END IF

!  if the string starts with a ., add a 0 at the front

      IF ( TRIM_VALUE( 1 : 1 ) == '.' ) THEN
        DO i = 24, 2, -1
          TRIM_VALUE( i : i ) = TRIM_VALUE( i - 1 : i - 1 )
        END DO
        TRIM_VALUE( 1 : 1 ) = '0'
      END IF

!  if the string starts with a -., replace by -0. at the front

      IF ( TRIM_VALUE( 1 : 2 ) == '-.' ) THEN
        DO i = 24, 3, -1
          TRIM_VALUE( i : i ) = TRIM_VALUE( i - 1 : i - 1 )
        END DO
        TRIM_VALUE( 2 : 2 ) = '0'
      END IF
      RETURN

!  end of function TRIM_VALUE

      END FUNCTION TRIM_VALUE

      FUNCTION TRIM_INT( i )
      CHARACTER ( LEN = 16 ) :: TRIM_INT
      INTEGER :: i

!  write integer as a left shifted length 16 character

      TRIM_INT = REPEAT( ' ', 16 )
      WRITE( TRIM_INT, "( I0 )" ) i
      RETURN

!  end of function TRIM_INT

      END FUNCTION TRIM_INT

!  data printing subroutines

      SUBROUTINE WRITE_X( out, n, X, X_l, X_u, Z )
      INTEGER :: n, out
      REAL ( KIND = wp ), DIMENSION( n ) :: X, X_l, X_u, Z
      WRITE( out, "( ' *       i      X_l          X          X_u          Z')")
      DO i = 1, n
        WRITE( out, "( ' * ', I7, 4ES12.4 )" )                                 &
          i, X_l( i ), X( i ), X_u( i ), Z( i )
      END DO
      END SUBROUTINE WRITE_X

      SUBROUTINE WRITE_Y( out, m, Y, C_l, C_u )
      INTEGER :: m, out
      REAL ( KIND = wp ), DIMENSION( m ) :: Y, C_l, C_u
      WRITE( out, "( ' *       i      C_l         C_u          Y   ' )" )
      DO i = 1, m
        WRITE( out, "( ' * ', I7, 3ES12.4 )" ) i, C_l( i ), C_u( i ), Y( i )
      END DO
      END SUBROUTINE WRITE_Y

      SUBROUTINE WRITE_X_type( out, n, X_type )
      INTEGER :: n, out
      INTEGER, DIMENSION( n ) :: X_type
      INTEGER :: i
      WRITE( out, "( ' *       i  X_type' )" )
      DO i = 1, n
        WRITE( out, "( ' * ', I7, 2X, I0 )" ) i, X_type( i )
      END DO
      END SUBROUTINE WRITE_X_type

      SUBROUTINE WRITE_p_name( out, p_name )
      INTEGER :: out
      CHARACTER ( len = 10 ) ::  p_name
      WRITE( out, "( ' * p_name = ', A )" ) p_name
      END SUBROUTINE WRITE_p_name

      SUBROUTINE WRITE_X_names( out, n, X_names )
      INTEGER :: n, out
      CHARACTER ( len = 10 ), DIMENSION( n ) :: X_names
      INTEGER :: i
      WRITE( out, "( ' *       i  X_name' )" )
      DO i = 1, n
        WRITE( out, "( ' * ', I7, 2X, A10 )" ) i, X_names( i )
      END DO
      END SUBROUTINE WRITE_X_names

      SUBROUTINE WRITE_C_names( out, m, C_names )
      INTEGER :: m, out
      CHARACTER ( len = 10 ), DIMENSION( m ) :: C_names
      INTEGER :: i
      WRITE( out, "( ' *       i  C_name' )" )
      DO i = 1, m
        WRITE( out, "( ' * ', I7, 2X, A10 )" ) i, C_names( i )
      END DO
      END SUBROUTINE WRITE_C_names

      SUBROUTINE WRITE_f( out, f )
      INTEGER :: out
      REAL ( KIND = wp ) :: f
      WRITE( out, "( ' * f = ', ES12.4 )" ) f
      END SUBROUTINE WRITE_f

      SUBROUTINE WRITE_C( out, m, C )
      INTEGER :: m, out
      REAL ( KIND = wp ), DIMENSION( m ) :: C
      INTEGER :: i
      WRITE( out, "( ' *       i       C' )" )
      DO i = 1, m
        WRITE( out, "( ' * ', I7, ES12.4 )" ) i, C( i )
      END DO
      END SUBROUTINE WRITE_C

      SUBROUTINE WRITE_G( out, n, G )
      INTEGER :: n, out
      REAL ( KIND = wp ), DIMENSION( n ) :: G
      INTEGER :: i
      WRITE( out, "( ' *       i       G' )" )
      DO i = 1, n
        WRITE( out, "( ' * ', I7, ES12.4 )" ) i, G( i )
      END DO
      END SUBROUTINE WRITE_G

      SUBROUTINE WRITE_H_dense( out, n, l_h2_1, H2_val )
      INTEGER :: n, l_h2_1, out
      REAL ( KIND = wp ), DIMENSION( l_h2_1, n ) :: H2_val
      INTEGER :: i, j
      WRITE( out, "( ' * H(dense)' )" )
      DO j = 1, n, 4
        IF ( j + 3 <= n ) THEN
          WRITE( out, "( ' *       i   j', I8, 3I12 )" ) j, j + 1, j + 2, j + 3
        ELSE IF ( j + 2 <= n ) THEN
          WRITE( out, "( ' *       i   j', I8, 2I12 )" ) j, j + 1, j + 2
        ELSE IF ( j + 1 <= n ) THEN
          WRITE( out, "( ' *       i   j', I8, I12 )" ) j, j + 1
        ELSE
          WRITE( out, "( ' *       i   j', I8 )" ) j
        END IF
        DO i = 1, n
          IF ( j + 3 <= n ) THEN
            WRITE( out, "( ' * ', I7,  4X, 4ES12.4 )" )                        &
              i, H2_val( i, j ), H2_val( i, j + 1 ),                           &
              H2_val( i, j + 2 ), H2_val( i, j + 3 )
          ELSE IF ( j + 2 <= n ) THEN
            WRITE( out, "( ' * ',  I7, 4X, 3ES12.4 )" )                        &
              i, H2_val( i, j ), H2_val( i, j + 1 ), H2_val( i, j + 2 )
          ELSE IF ( j + 1 <= n ) THEN
            WRITE( out, "( ' * ',  I7, 4X, 2ES12.4 )" )                        &
              i, H2_val( i, j ), H2_val( i, j + 1 )
          ELSE
            WRITE( out, "( ' * ',  I7, 4X, ES12.4 )" ) i, H2_val( i, j )
          END IF
        END DO
      END DO
      END SUBROUTINE WRITE_H_dense

      SUBROUTINE WRITE_A_dense( out, n, m, l_j2_1, l_j2_2, J2_val )
      INTEGER :: n, m, l_J2_1, out
      REAL ( KIND = wp ), DIMENSION( l_j2_1, l_j2_2 ) :: J2_val
      INTEGER :: i, j
      WRITE( out, "( ' * A(dense)' )" )
      DO j = 1, n, 4
        IF ( j + 3 <= n ) THEN
          WRITE( out, "( ' *       i   j', I8, 3I12 )" ) j, j + 1, j + 2, j + 3
        ELSE IF ( j + 2 <= n ) THEN
          WRITE( out, "( ' *       i   j', I8, 2I12 )" ) j, j + 1, j + 2
        ELSE IF ( j + 1 <= n ) THEN
          WRITE( out, "( ' *       i   j', I8, I12 )" ) j, j + 1
        ELSE
          WRITE( out, "( ' *       i   j', I8 )" ) j
        END IF
        DO i = 1, m
          IF ( j + 3 <= n ) THEN
            WRITE( out, "( ' * ', I7,  4X, 4ES12.4 )" )                        &
              i, J2_val( i, j ), J2_val( i, j + 1 ),                           &
              J2_val( i, j + 2 ), J2_val( i, j + 3 )
          ELSE IF ( j + 2 <= n ) THEN
            WRITE( out, "( ' * ',  I7, 4X, 3ES12.4 )" )                        &
              i, J2_val( i, j ), J2_val( i, j + 1 ), J2_val( i, j + 2 )
          ELSE IF ( j + 1 <= n ) THEN
            WRITE( out, "( ' * ',  I7, 4X, 2ES12.4 )" )                        &
              i, J2_val( i, j ), J2_val( i, j + 1 )
          ELSE
            WRITE( out, "( ' * ',  I7, 4X, ES12.4 )" ) i, J2_val( i, j )
          END IF
        END DO
      END DO
      END SUBROUTINE WRITE_A_dense

      SUBROUTINE WRITE_H_sparse( out, H_ne, l_h, H_val, H_row, H_col )
      INTEGER :: l_h, H_ne, out
      INTEGER, DIMENSION( l_h ) :: H_row, H_col
      REAL ( KIND = wp ), DIMENSION( l_h ) :: H_val
      INTEGER :: i
      IF ( H_ne == 0 ) RETURN
      WRITE( out, "( ' * H(sparse)' )" )
      WRITE( out, "( ' * ', 2( '    row    col     val    ' ) )" )
      DO i = 1, H_ne, 2
        IF ( i + 1 <= H_ne ) THEN
          WRITE( out, "( ' * ',  2( 2I7, ES12.4 ) )" )                         &
            H_row( i ), H_col( i ), H_val( i ),                                &
            H_row( i + 1 ), H_col( i + 1 ), H_val( i + 1 )
        ELSE
          WRITE( out, "( ' * ',  2( 2I7, ES12.4 ) )" )                         &
            H_row( i ), H_col( i ), H_val( i )
        END IF
      END DO
      END SUBROUTINE WRITE_H_sparse

      SUBROUTINE WRITE_A_sparse( out, A_ne, l_a, A_val, A_row, A_col )
      INTEGER :: l_a, A_ne, out
      INTEGER, DIMENSION( l_a ) :: A_row, A_col
      REAL ( KIND = wp ), DIMENSION( l_a ) :: A_val
      INTEGER :: i
      IF ( A_ne == 0 ) RETURN
      WRITE( out, "( ' * A(sparse)' )" )
      WRITE( out, "( ' * ', 2( '    row    col     val    ' ) )" )
      DO i = 1, A_ne, 2
        IF ( i + 1 <= A_ne ) THEN
          WRITE( out, "( ' * ',  2( 2I7, ES12.4 ) )" )                         &
            A_row( i ), A_col( i ), A_val( i ),                                &
            A_row( i + 1 ), A_col( i + 1 ), A_val( i + 1 )
        ELSE
          WRITE( out, "( ' * ',  2( 2I7, ES12.4 ) )" )                         &
            A_row( i ), A_col( i ), A_val( i )
        END IF
      END DO
      END SUBROUTINE WRITE_A_sparse

      SUBROUTINE WRITE_H_byrows( out, n, H_val, H_col, H_ptr )
      INTEGER :: n, out
      INTEGER, DIMENSION( n + 1 ) :: H_ptr
      INTEGER, DIMENSION( H_ptr( n + 1 ) - 1 ) :: H_col
      REAL ( KIND = wp ), DIMENSION( H_ptr( n + 1 ) - 1 ) :: H_val
      INTEGER :: i, l_up, maxc
      WRITE( out, "( ' * H(by rows)' )" )
      maxc = MAXVAL( H_ptr( 2 : n + 1 ) - H_ptr( 1 : n ) )
      IF ( maxc >= 3 ) THEN
        WRITE( out, "( ' *     row ', 3( '   col     val     ' ) )" )
      ELSE IF ( maxc >= 2 ) THEN
        WRITE( out, "( ' *     row ', 2( '   col     val     ' ) )" )
      ELSE IF ( maxc >= 1 ) THEN
        WRITE( out, "( ' *     row ',  ( '   col     val     ' ) )" )
      ELSE
        WRITE( out, "( ' *     row ' )" )
      END IF
      DO i = 1, n
        l_up =  H_ptr( i + 1 ) - 1
        IF ( H_ptr( i ) > l_up ) THEN
          WRITE( out, "( ' * ',  I7, ' no entries ' )" ) i
        ELSE
          DO l = H_ptr( i ), l_up, 3
            IF ( l + 2 <= l_up ) THEN
              WRITE( out, "( ' * ',  I7, 3( I7, ES12.4 ) )" ) i,               &
                H_col( l ), H_val( l ), H_col( l + 1 ), H_val( l + 1 ),        &
                H_col( l + 2 ), H_val( l + 2 )
            ELSE IF ( l + 1 <= l_up ) THEN
              WRITE( out, "( ' * ',  I7, 2( I7, ES12.4 ) )" ) i,               &
                H_col( l ), H_val( l ), H_col( l + 1 ), H_val( l + 1 )
            ELSE
              WRITE( out, "( ' * ',  I7, ( I7, ES12.4 ) )" ) i,                &
                H_col( l ), H_val( l )
            END IF
          END DO
        END IF
      END DO
      END SUBROUTINE WRITE_H_byrows

      SUBROUTINE WRITE_A_byrows( out, m, A_val, A_col, A_ptr )
      INTEGER :: m, out
      INTEGER, DIMENSION( m + 1 ) :: A_ptr
      INTEGER, DIMENSION( A_ptr( m + 1 ) - 1 ) :: A_col
      REAL ( KIND = wp ), DIMENSION( A_ptr( m + 1 ) - 1 ) :: A_val
      INTEGER :: i, l_up, maxc
      WRITE( out, "( ' * A(by rows)' )" )
      maxc = MAXVAL( A_ptr( 2 : m + 1 ) - A_ptr( 1 : m ) )
      IF ( maxc >= 3 ) THEN
        WRITE( out, "( ' *     row ', 3( '   col     val     ' ) )" )
      ELSE IF ( maxc >= 2 ) THEN
        WRITE( out, "( ' *     row ', 2( '   col     val     ' ) )" )
      ELSE IF ( maxc >= 1 ) THEN
        WRITE( out, "( ' *     row ',  ( '   col     val     ' ) )" )
      ELSE
        WRITE( out, "( ' *     row ' )" )
      END IF
      DO i = 1, m
        l_up =  A_ptr( i + 1 ) - 1
        IF ( A_ptr( i ) > l_up ) THEN
          WRITE( out, "( ' * ',  I7, ' no entries ' )" ) i
        ELSE
          DO l = A_ptr( i ), l_up, 3
            IF ( l + 2 <= l_up ) THEN
              WRITE( out, "( ' * ',  I7, 3( I7, ES12.4 ) )" ) i,               &
                A_col( l ), A_val( l ), A_col( l + 1 ), A_val( l + 1 ),        &
                A_col( l + 2 ), A_val( l + 2 )
            ELSE IF ( l + 1 <= l_up ) THEN
              WRITE( out, "( ' * ',  I7, 2( I7, ES12.4 ) )" ) i,               &
                A_col( l ), A_val( l ), A_col( l + 1 ), A_val( l + 1 )
            ELSE
              WRITE( out, "( ' * ',  I7, ( I7, ES12.4 ) )" ) i,                &
                A_col( l ), A_val( l )
            END IF
          END DO
        END IF
      END DO
      END SUBROUTINE WRITE_A_byrows

      SUBROUTINE WRITE_H_bycols( out, n, H_val, H_row, H_ptr )
      INTEGER :: n, out
      INTEGER, DIMENSION( n + 1 ) :: H_ptr
      INTEGER, DIMENSION( H_ptr( n + 1 ) - 1 ) :: H_row
      REAL ( KIND = wp ), DIMENSION( H_ptr( n + 1 ) - 1 ) :: H_val
      INTEGER :: i, l_up, maxr
      WRITE( out, "( ' * H(by cols)' )" )
      maxr = MAXVAL( H_ptr( 2 : n + 1 ) - H_ptr( 1 : n ) )
      IF ( maxr >= 3 ) THEN
        WRITE( out, "( ' *     col ', 3( '   row     val     ' ) )" )
      ELSE IF ( maxr >= 2 ) THEN
        WRITE( out, "( ' *     col ', 2( '   row     val     ' ) )" )
      ELSE IF ( maxr >= 1 ) THEN
        WRITE( out, "( ' *     col ',  ( '   row     val     ' ) )" )
      ELSE
        WRITE( out, "( ' *     col ' )" )
      END IF
      DO i = 1, n
        l_up =  H_ptr( i + 1 ) - 1
        IF ( H_ptr( i ) > l_up ) THEN
          WRITE( out, "( ' * ',  I7, ' no entries ' )" ) i
        ELSE
          DO l = H_ptr( i ), l_up, 3
            IF ( l + 2 <= l_up ) THEN
              WRITE( out, "( ' * ',  I7, 3( I7, ES12.4 ) )" ) i,               &
                H_row( l ), H_val( l ), H_row( l + 1 ), H_val( l + 1 ),        &
                H_row( l + 2 ), H_val( l + 2 )
            ELSE IF ( l + 1 <= l_up ) THEN
              WRITE( out, "( ' * ',  I7, 2( I7, ES12.4 ) )" ) i,               &
                H_row( l ), H_val( l ), H_row( l + 1 ), H_val( l + 1 )
            ELSE
              WRITE( out, "( ' * ',  I7, ( I7, ES12.4 ) )" ) i,                &
                H_row( l ), H_val( l )
            END IF
          END DO
        END IF
      END DO
      END SUBROUTINE WRITE_H_bycols

      SUBROUTINE WRITE_A_bycols( out, n, A_val, A_row, A_ptr )
      INTEGER :: n, out
      INTEGER, DIMENSION( n + 1 ) :: A_ptr
      INTEGER, DIMENSION( A_ptr( n + 1 ) - 1 ) :: A_row
      REAL ( KIND = wp ), DIMENSION( A_ptr( n + 1 ) - 1 ) :: A_val
      INTEGER :: i, l_up, maxr
      WRITE( out, "( ' * A(by cols)' )" )
      maxr = MAXVAL( A_ptr( 2 : n + 1 ) - A_ptr( 1 : n ) )
      IF ( maxr >= 3 ) THEN
        WRITE( out, "( ' *     col ', 3( '   row     val     ' ) )" )
      ELSE IF ( maxr >= 2 ) THEN
        WRITE( out, "( ' *     col ', 2( '   row     val     ' ) )" )
      ELSE IF ( maxr >= 1 ) THEN
        WRITE( out, "( ' *     col ',  ( '   row     val     ' ) )" )
      ELSE
        WRITE( out, "( ' *     col ' )" )
      END IF
      DO i = 1, n
        l_up =  A_ptr( i + 1 ) - 1
        IF ( A_ptr( i ) > l_up ) THEN
          WRITE( out, "( ' * ',  I7, ' no entries ' )" ) i
        ELSE
          DO l = A_ptr( i ), l_up, 3
            IF ( l + 2 <= l_up ) THEN
              WRITE( out, "( ' * ',  I7, 3( I7, ES12.4 ) )" ) i,               &
                A_row( l ), A_val( l ), A_row( l + 1 ), A_val( l + 1 ),        &
                A_row( l + 2 ), A_val( l + 2 )
            ELSE IF ( l + 1 <= l_up ) THEN
              WRITE( out, "( ' * ',  I7, 2( I7, ES12.4 ) )" ) i,               &
                A_row( l ), A_val( l ), A_row( l + 1 ), A_val( l + 1 )
            ELSE
              WRITE( out, "( ' * ',  I7, ( I7, ES12.4 ) )" ) i,                &
                A_row( l ), A_val( l )
            END IF
          END DO
        END IF
      END DO
      END SUBROUTINE WRITE_A_bycols

    END PROGRAM CUTEST_qplib_main

