! THIS VERSION: CUTEST 1.0 - 17/11/2012 AT 10:00 GMT.

!- C U T E S T  t e s t _ u n c o n s t r a i n e d _ t o o l s  P R O G R A M -

    PROGRAM CUTEST_test_unconstrained_tools

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 2003 version released November 2012

      USE CUTEST

!--------------------
!   P r e c i s i o n
!--------------------

      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!----------------------
!   P a r a m e t e r s
!----------------------

      INTEGER, PARAMETER :: input = 55
      INTEGER, PARAMETER :: out = 6
      INTEGER, PARAMETER :: buffer = 77
      REAL ( KIND = wp ), PARAMETER :: one = 1.0_wp

!--------------------------------
!   L o c a l   V a r i a b l e s
!--------------------------------

      INTEGER :: n, HE_nel, HE_val_ne, HE_row_ne, status, alloc_stat
      INTEGER :: l_h2_1, l_h, l_hel, H_ne, l_he_val, l_he_row
      REAL ( KIND = wp ) :: f
      LOGICAL :: grad, byrows, goth
      CHARACTER ( len = 10 ) ::  p_name
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: X_type, H_row, H_col
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: HE_row, HE_row_ptr, HE_val_ptr
      REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: X, X_l, X_u, G
      REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: H_val, HE_val, P, HP
      REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : , : ) :: H2_val, H_band
      CHARACTER ( len = 10 ), ALLOCATABLE, DIMENSION( : ) :: X_names
!     REAL :: CPU( 2 ), CALLS( 4 )
      TYPE ( CUTEST_data_type ) :: data

!  open the problem data file

      OPEN ( input, FILE = 'u_OUTSDIF.d', FORM = 'FORMATTED', STATUS = 'OLD' )

!  allocate basic arrays

      WRITE( out, "( ' Call UDIMEN ' )" )
      CALL UDIMEN( input, status, n )
      l_h2_1 = n
      ALLOCATE( X( n ), X_l( n ), X_u( n ), G( n ), P( n ), HP( n ),           &
                X_names( n ), X_type( n ), stat = alloc_stat )
      IF ( alloc_stat /= 0 ) GO TO 990
      ALLOCATE( H2_val( l_h2_1, n ), stat = alloc_stat )
      IF ( alloc_stat /= 0 ) GO TO 990

!  set up SIF data

      WRITE( out, "( ' Call USETUP ' )" )
      CALL USETUP( data, status, input, out, buffer, n, X, X_l, X_u )
      IF ( status /= 0 ) GO to 900

!  obtain variable and problem names

      WRITE( out, "( ' Call UNAMES' )" )
      CALL UNAMES( data, status, n, p_name, X_names )
      IF ( status /= 0 ) GO to 900

!  obtain problem name

      WRITE( out, "( ' Call PBNAME' )" )
      CALL PBNAME( data, status, p_name )
      IF ( status /= 0 ) GO to 900

!  obtain variable names

      WRITE( out, "( ' Call VARNAMES' )" )
      CALL VARNAMES( data, status, n, X_names )
      IF ( status /= 0 ) GO to 900

!  obtain variable types

      WRITE( out, "( ' Call UVARTY' )" )
      CALL UVARTY( data, status, n, X_type )
      IF ( status /= 0 ) GO to 900

!  compute the objective function value

      WRITE( out, "( ' Call UFN' )" )
      CALL UFN( data, status, n, X, f )
      IF ( status /= 0 ) GO to 900

!  compute the gradient value

      WRITE( out, "( ' Call UGR' )" )
      CALL UGR( data, status, n, X, G )
      IF ( status /= 0 ) GO to 900

!  compute the objective function and gradient values

      grad = .TRUE.
      WRITE( out, "( ' Call UOFG with grad = .TRUE.' )" )
      CALL UOFG( data, status, n, X, f, G, grad )
      IF ( status /= 0 ) GO to 900

      grad = .FALSE.
      WRITE( out, "( ' Call UOFG with grad = .FALSE.' )" )
      CALL UOFG( data, status, n, X, f, G, grad )
      IF ( status /= 0 ) GO to 900

!  compute the dense Hessian value

      WRITE( out, "( ' Call UDH' )" )
      CALL UDH( data, status, n, X, l_h2_1, H2_val )
      IF ( status /= 0 ) GO to 900

!  compute the gradient and dense Hessian values

      WRITE( out, "( ' Call UGRDH' )" )
      CALL UGRDH( data, status, n, X, G, l_h2_1, H2_val )
      IF ( status /= 0 ) GO to 900

!  compute the number of nonzeros in the sparse Hessian

      WRITE( out, "( ' Call UDIMSH' )" )
      CALL UDIMSH( data, status, H_ne )
      IF ( status /= 0 ) GO to 900

      l_h = H_ne
      ALLOCATE( H_val( l_h ), H_row( l_h ), H_col( l_h ), stat = alloc_stat )
      IF ( alloc_stat /= 0 ) GO TO 990

!  compute the sparse Hessian value

      WRITE( out, "( ' Call USH' )" )
      CALL USH( data, status, n, X, H_ne, l_h, H_val, H_row, H_col )
      IF ( status /= 0 ) GO to 900

!  compute the gradient and sparse Hessian values

      WRITE( out, "( ' Call UGRSH' )" )
      CALL UGRSH( data, status, n, X, G, H_ne, l_h, H_val, H_row, H_col )
      IF ( status /= 0 ) GO to 900

!  compute the number of nonzeros in the element Hessian

      WRITE( out, "( ' Call UDIMSE' )" )
      CALL UDIMSE( data, status, HE_nel, HE_val_ne, HE_row_ne )
      IF ( status /= 0 ) GO to 900

      l_hel = HE_nel + 1
      l_he_val = HE_val_ne
      l_he_row = HE_row_ne
      ALLOCATE( HE_row_ptr( l_hel ), HE_val_ptr( l_hel ), HE_row( l_he_row ),  &
                HE_val( l_he_val ), stat = alloc_stat )
      IF ( alloc_stat /= 0 ) GO TO 990

!  compute the element Hessian value

      byrows = .FALSE.
      WRITE( out, "( ' Call UEH with byrows = .FALSE.' )" )
      CALL UEH( data, status, n, X, HE_nel, HE_row, l_he_row, l_hel,           &
                HE_row_ptr, HE_val, l_he_val, HE_val_ptr, byrows )
      IF ( status /= 0 ) GO to 900
      byrows = .TRUE.
      WRITE( out, "( ' Call UEH with byrows = .TRUE.' )" )
      CALL UEH( data, status, n, X, HE_nel, HE_row, l_he_row, l_hel,           &
                HE_row_ptr, HE_val, l_he_val, HE_val_ptr, byrows )
      IF ( status /= 0 ) GO to 900

!  compute the gradient and element Hessian values

      byrows = .FALSE.
      WRITE( out, "( ' Call UGREH with byrows = .FALSE' )" )
      CALL UGREH( data, status, n, X, G, HE_nel, HE_row, l_he_row, l_hel,       &
                  HE_row_ptr, HE_val, l_he_val, HE_val_ptr, byrows )
      IF ( status /= 0 ) GO to 900
      byrows = .TRUE.
      WRITE( out, "( ' Call UGREH with byrows = .TRUE.' )" )
      CALL UGREH( data, status, n, X, G, HE_nel, HE_row, l_he_row, l_hel,       &
                  HE_row_ptr, HE_val, l_he_val, HE_val_ptr, byrows )
      IF ( status /= 0 ) GO to 900

!  compute a Hessian-vector product

      P = one
      goth = .FALSE.
      WRITE( out, "( ' Call UPROD with goth = .FALSE.' )" )
      CALL UPROD( data, status, n, goth, X, P, HP )
      IF ( status /= 0 ) GO to 900
      goth = .TRUE.
      WRITE( out, "( ' Call UPROD with goth = .TRUE.' )" )
      CALL UPROD( data, status, n, goth, X, P, HP )
      IF ( status /= 0 ) GO to 900

!  compute a band of the Hessian

      nsemib = n / 2
      lbandh = nsemib
      ALLOCATE( H_band( 0 : lbandh, n ), stat = alloc_stat )
      IF ( alloc_stat /= 0 ) GO TO 990

      goth = .FALSE.
      WRITE( out, "( ' Call UBANDH with goth = .FALSE.' )" )
      CALL UBANDH( data, status, n, goth, X, nsemib, H_band, lbandh )
      IF ( status /= 0 ) GO to 900
      goth = .TRUE.
      WRITE( out, "( ' Call UBANDH with goth = .TRUE.' )" )
      CALL UBANDH( data, status, n, goth, X, nsemib, H_band, lbandh )
      IF ( status /= 0 ) GO to 900

!  terminal exit

!     WRITE( out, "( ' Call UREPRT' )" )
!     CALL UREPRT( CALLS, CPU )

      DEALLOCATE( X_type, H_row, H_col, HE_row, HE_row_ptr, HE_val_ptr, X,     &
                  X_l, X_u, G, H_val, HE_val, P, HP, H2_val, H_band, X_names,  &
                  stat = alloc_stat )
      CLOSE( input )
      STOP

!  error exits

 900  CONTINUE
      WRITE( out, "( ' error status = ', I0 )" ) status
      CLOSE( INPUT  )
      STOP

 990  CONTINUE
      WRITE( out, "( ' Allocation error, status = ', I0 )" ) alloc_stat
      CLOSE( INPUT  )
      STOP

    END PROGRAM CUTEST_test_unconstrained_tools
