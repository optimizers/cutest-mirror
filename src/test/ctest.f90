! THIS VERSION: CUTEST 1.0 - 17/11/2012 AT 17:00 GMT.

!-*- C U T E S T  t e s t _ c o n s t r a i n e d _ t o o l s  P R O G R A M -*-

    PROGRAM CUTEST_test_constrained_tools

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
      REAL ( KIND = wp ), PARAMETER :: one = 1.0_wp

!--------------------------------
!   L o c a l   V a r i a b l e s
!--------------------------------

      INTEGER :: n, m, H_ne, HE_nel, HE_val_ne, HE_row_ne, J_ne, Ji_ne, status
      INTEGER :: l_h2_1, l_h, l_hel, l_he_val, l_he_row, alloc_stat
      REAL ( KIND = wp ) :: f, ci
      LOGICAL :: grad, byrows, goth, efirst, lfirst, nvfrst, grlagf, jtrans
      CHARACTER ( len = 10 ) ::  p_name
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: H_row, H_col, X_type
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: HE_row, HE_row_ptr, HE_val_ptr
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: J_row, J_col
      REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: X, X_l, X_u, G, Ji
      REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: Y, C_l, C_u, C, J_val
      REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: H_val, HE_val, P, HP
      REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : , : ) :: H2_val, H_band
      REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : , : ) :: J2_val
      LOGICAL, ALLOCATABLE, DIMENSION( : ) :: EQUATION, LINEAR
      CHARACTER ( len = 10 ), ALLOCATABLE, DIMENSION( : ) :: X_names, C_names
!     REAL :: CPU( 2 ), CALLS( 4 )
      TYPE ( CUTEST_data_type ) :: data

!  open the problem data file

      OPEN ( input, FILE = 'c_OUTSDIF.d', FORM = 'FORMATTED', STATUS = 'OLD' )

!  allocate basic arrays

      WRITE( out, "( ' CALL CDIMEN ' )" )
      CALL CDIMEN( input, status, n, m )
      IF ( alloc_stat /= 0 ) GO TO 990
      l_h2_1 = n
      ALLOCATE( X( n ), X_l( n ), X_u( n ), G( n ), P( n ), Ji( n ),           &
                HP( n ), X_names( n ), X_type( n ), stat = alloc_stat )
      IF ( alloc_stat /= 0 ) GO TO 990
      ALLOCATE( C( m ), Y( m ), C_l( m ), C_u( m ), C_names( m ),              &
                EQUATION( m ), LINEAR( m ), stat = alloc_stat )
      IF ( alloc_stat /= 0 ) GO TO 990
      ALLOCATE( H2_val( l_h2_1, n ), stat = alloc_stat )
      IF ( alloc_stat /= 0 ) GO TO 990
      l_j2_1 = MAX( m, n ) ; l_j2_2 = l_j2_1      
      ALLOCATE( J2_val( l_j2_1, l_j2_2 ), stat = alloc_stat )
      IF ( alloc_stat /= 0 ) GO TO 990

!  set up SIF data

      efirst = .TRUE. ; lfirst = .TRUE. ; nvfrst = .TRUE.
      WRITE( out, "( ' CALL CSETUP ' )" )
      CALL CSETUP( data, status, input, out, n, m, X, X_l, X_u,                 &
                   EQUATION, LINEAR, Y, C_l, C_u, efirst, lfirst, nvfrst )
      IF ( status /= 0 ) GO to 900

!  obtain variable and problem names

      WRITE( out, "( ' CALL CNAMES' )" )
      CALL CNAMES( data, status, n, m, p_name, X_names, C_names )
      IF ( status /= 0 ) GO to 900

!  obtain constraint names

      WRITE( out, "( ' Call CONNAMES' )" )
      CALL CONNAMES( data, status, m, C_names )
      IF ( status /= 0 ) GO to 900

!  obtain variable types

      WRITE( out, "( ' CALL CVARTY' )" )
      CALL CVARTY( data, status, n, X_type )
      IF ( status /= 0 ) GO to 900

!  compute the objective function value

      WRITE( out, "( ' CALL CFN' )" )
      CALL CFN( data, status, n, m, X, f, C )
      IF ( status /= 0 ) GO to 900

!  compute the gradient and dense Jacobian values

      grlagf = .TRUE. ; jtrans = .TRUE.
      WRITE( out, "( ' CALL CGR with grlagf = .TRUE. and jtrans = .TRUE.' )" )
      CALL CGR( data, status, n, m, X, grlagf, m, Y, G, jtrans,                 &
                l_j2_1, l_j2_2, J2_val )
      IF ( status /= 0 ) GO to 900
      grlagf = .TRUE. ; jtrans = .FALSE.
      WRITE( out, "( ' CALL CGR with grlagf = .TRUE. and jtrans = .FALSE.' )" )
      CALL CGR( data, status, n, m, X, grlagf, m, Y, G, jtrans,                 &
                l_j2_1, l_j2_2, J2_val )
      IF ( status /= 0 ) GO to 900
      grlagf = .FALSE. ; jtrans = .TRUE.
      WRITE( out, "( ' CALL CGR with grlagf = .FALSE. and jtrans = .TRUE.' )" )
      CALL CGR( data, status, n, m, X, grlagf, m, Y, G, jtrans,                 &
                l_j2_1, l_j2_2, J2_val )
      IF ( status /= 0 ) GO to 900
      grlagf = .FALSE. ; jtrans = .FALSE.
      WRITE( out, "( ' CALL CGR with grlagf = .FALSE. and jtrans = .FALSE.' )" )
      CALL CGR( data, status, n, m, X, grlagf, m, Y, G, jtrans,                 &
                l_j2_1, l_j2_2, J2_val )
      IF ( status /= 0 ) GO to 900

!  compute the objective function and gradient values

      grad = .TRUE.
      WRITE( out, "( ' CALL COFG with grad = .TRUE.' )" )
      CALL COFG( data, status, n, X, f, G, grad )
      IF ( status /= 0 ) GO to 900
      grad = .FALSE.
      WRITE( out, "( ' CALL COFG with grad = .FALSE.' )" )
      CALL COFG( data, status, n, X, f, G, grad )
      IF ( status /= 0 ) GO to 900

!  compute the number of nonzeros in the sparse Jacobian

      WRITE( out, "( ' CALL CDIMSJ' )" )
      CALL CDIMSJ( data, status, J_ne )
      IF ( status /= 0 ) GO to 900

      l_j = J_ne
      ALLOCATE( J_val( l_j ), J_row( l_j ), J_col( l_j ), stat = alloc_stat )
      IF ( alloc_stat /= 0 ) GO TO 990

!  compute the gradient and sparse Jacobian values

      grlagf = .TRUE.
      WRITE( out, "( ' CALL CSGR with grlagf = .TRUE.' )" )
      CALL CSGR( data, status, n, m, grlagf, m, Y, X, J_ne, l_j,               &
                 J_val, J_col, J_row )
      IF ( status /= 0 ) GO to 900
      grlagf = .FALSE.
      WRITE( out, "( ' CALL CSGR with grlagf = .FALSE.' )" )
      CALL CSGR( data, status, n, m, grlagf, m, Y, X, J_ne, l_j,               &
                 J_val, J_col, J_row )
      IF ( status /= 0 ) GO to 900

!  compute the constraint and dense Jacobian values

      grad = .TRUE. ; jtrans = .TRUE.
      WRITE( out, "( ' CALL CCFG with grad = .TRUE. and jtrans = .TRUE.' )" )
      CALL CCFG( data, status, n, m, X, C, jtrans, l_j2_1, l_j2_2, J_val, grad )
      IF ( status /= 0 ) GO to 900
      grad = .TRUE. ; jtrans = .FALSE.
      WRITE( out, "( ' CALL CCFG with grad = .TRUE. and jtrans = .FALSE.' )" )
      CALL CCFG( data, status, n, m, X, C, jtrans, l_j2_1, l_j2_2, J_val, grad )
      IF ( status /= 0 ) GO to 900
      grad = .FALSE. ; jtrans = .TRUE.
      WRITE( out, "( ' CALL CCFG with grad = .FALSE. and jtrans = .TRUE.' )" )
      CALL CCFG( data, status, n, m, X, C, jtrans, l_j2_1, l_j2_2, J_val, grad )
      IF ( status /= 0 ) GO to 900
      grad = .FALSE. ; jtrans = .FALSE.
      WRITE( out, "( ' CALL CCFG with grad = .FALSE. and jtrans = .FALSE.' )" )
      CALL CCFG( data, status, n, m, X, C, jtrans, l_j2_1, l_j2_2, J_val, grad )
      IF ( status /= 0 ) GO to 900

!  compute the constraint and sparse Jacobian values

      grad = .TRUE.
      WRITE( out, "( ' CALL CSCFG with grad = .TRUE.' )" )
      CALL CSCFG( data, status, n, m, X, m, C, J_ne, l_j, J_val, J_col, J_row, &
                  grad )
      IF ( status /= 0 ) GO to 900
      grad = .FALSE.
      WRITE( out, "( ' CALL CSCFG with grad = .FALSE.' )" )
      CALL CSCFG( data, status, n, m, X, m, C, J_ne, l_j, J_val, J_col, J_row, &
                  grad )
      IF ( status /= 0 ) GO to 900

!  compute an individual constraint and its dense gradient

      icon = 1
      grad = .TRUE.
      WRITE( out, "( ' CALL CCIFG with grad = .TRUE.' )" )
      CALL CCIFG( data, status, n, icon, X, ci, Ji, grad )
      IF ( status /= 0 ) GO to 900
      grad = .FALSE.
      WRITE( out, "( ' CALL CCIFG with grad = .FALSE.' )" )
      CALL CCIFG( data, status, n, icon, X, ci, Ji, grad )
      IF ( status /= 0 ) GO to 900

!  compute an individual constraint and its dense gradient

      grad = .TRUE.
      WRITE( out, "( ' CALL CSCIFG with grad = .TRUE.' )" )
      CALL CSCIFG( data, status, n, icon, X, ci, Ji_ne, n, Ji, J_col, grad )
      IF ( status /= 0 ) GO to 900
      grad = .FALSE.
      WRITE( out, "( ' CALL CSCIFG with grad = .FALSE.' )" )
      CALL CSCIFG( data, status, n, icon, X, ci, Ji_ne, n, Ji, J_col, grad )
      IF ( status /= 0 ) GO to 900

!  compute the dense Hessian value

      WRITE( out, "( ' CALL CDH' )" )
      CALL CDH( data, status, n, m, X, m, Y, l_h2_1, H2_val )
      IF ( status /= 0 ) GO to 900

!  compute the dense Hessian value of the objective or a constraint

      iprob = 0
      WRITE( out, "( ' CALL CIDH for objective' )" )
      CALL CIDH( data, status, n, X, iprob, l_h2_1, H2_val )
      IF ( status /= 0 ) GO to 900
      iprob = 1
      WRITE( out, "( ' CALL CIDH for a constraint' )" )
      CALL CIDH( data, status, n, X, iprob, l_h2_1, H2_val )
      IF ( status /= 0 ) GO to 900

!  compute the gradient and dense Hessian values

      grlagf = .TRUE. ; jtrans = .TRUE.
      WRITE( out, "( ' CALL CGRDH with grlagf = .TRUE. and jtrans = .TRUE.' )" )
      CALL CGRDH( data, status, n, m, X, grlagf, m, Y, G, jtrans,              &
                  l_c2_1, l_c2_2, J2_val, l_h2_1, H2_val )
      IF ( status /= 0 ) GO to 900
      grlagf = .TRUE. ; jtrans = .FALSE.
      WRITE( out, "( ' CALL CGRDH with grlagf = .TRUE. and jtrans = .FALSE.' )")
      CALL CGRDH( data, status, n, m, X, grlagf, m, Y, G, jtrans,              &
                  l_c2_1, l_c2_2, J2_val, l_h2_1, H2_val )
      IF ( status /= 0 ) GO to 900
      grlagf = .FALSE. ; jtrans = .TRUE.
      WRITE( out, "( ' CALL CGRDH with grlagf = .FALSE. and jtrans = .TRUE.' )")
      CALL CGRDH( data, status, n, m, X, grlagf, m, Y, G, jtrans,              &
                  l_c2_1, l_c2_2, J2_val, l_h2_1, H2_val )
      IF ( status /= 0 ) GO to 900
      grlagf = .FALSE. ; jtrans = .FALSE.
      WRITE( out, "( ' CALL CGRDH with grlagf = .FALSE. and jtrans = .FALSE.')")
      CALL CGRDH( data, status, n, m, X, grlagf, m, Y, G, jtrans,              &
                  l_c2_1, l_c2_2, J2_val, l_h2_1, H2_val )
      IF ( status /= 0 ) GO to 900

!  compute the number of nonzeros in the sparse Hessian

      WRITE( out, "( ' CALL CDIMSH' )" )
      CALL CDIMSH( data, status, H_ne )
      IF ( status /= 0 ) GO to 900

      l_h = H_ne
      ALLOCATE( H_val( l_h ), H_row( l_h ), H_col( l_h ), stat = alloc_stat )
      IF ( alloc_stat /= 0 ) GO TO 990

!  compute the sparse Hessian value

      WRITE( out, "( ' CALL CSH' )" )
      CALL CSH( data, status, n, m, X, m, Y, H_ne, l_h, H_val, H_row, H_col )
      IF ( status /= 0 ) GO to 900

!  compute the sparse Hessian value without the objective

      WRITE( out, "( ' CALL CSH1' )" )
      CALL CSH1( data, status, n, m, X, m, Y, H_ne, l_h, H_val, H_row, H_col )
      IF ( status /= 0 ) GO to 900

!  compute the sparse Hessian value of the objective or a constraint

      iprob = 0
      WRITE( out, "( ' CALL CIDH for objective' )" )
      CALL CISH( data, status, n, X, iprob, H_ne, l_h, H_val, H_row, H_col )
      IF ( status /= 0 ) GO to 900
      iprob = 1
      WRITE( out, "( ' CALL CIDH for a constraint' )" )
      CALL CISH( data, status, n, X, iprob, H_ne, l_h, H_val, H_row, H_col )
      IF ( status /= 0 ) GO to 900

!  compute the gradient and sparse Hessian values

      grlagf = .TRUE.
      WRITE( out, "( ' CALL CSGRSH with grlagf = .TRUE.' )" )
      CALL CSGRSH( data, status, n, m, X, grlagf, m, Y, J_ne, l_j, J_val,      &
                   J_col, J_row, H_ne, l_h, H_val, H_row, H_col )
      IF ( status /= 0 ) GO to 900
      grlagf = .FALSE.
      WRITE( out, "( ' CALL CSGRSH with grlagf = .FALSE.' )" )
      CALL CSGRSH( data, status, n, m, X, grlagf, m, Y, J_ne, l_j, J_val, 
                   J_col, J_row, H_ne, l_h, H_val, H_row, H_col )
      IF ( status /= 0 ) GO to 900

!  compute the number of nonzeros in the element Hessian

      WRITE( out, "( ' CALL CDIMSE' )" )
      CALL CDIMSE( data, status, HE_nel, HE_val_ne, HE_row_ne )
      IF ( status /= 0 ) GO to 900

      l_hel = HE_nel + 1
      l_he_val = HE_val_ne
      l_he_row = HE_row_ne
      ALLOCATE( HE_row_ptr( l_hel ), HE_val_ptr( l_hel ), HE_row( l_he_row ),  &
                HE_val( l_he_val ), stat = alloc_stat )
      IF ( alloc_stat /= 0 ) GO TO 990

!  compute the element Hessian value

      byrows = .FALSE.
      WRITE( out, "( ' CALL CEH with byrows = .FALSE.' )" )
      CALL CEH( data, status, n, m, X, m, Y, HE_nel, HE_row, l_he_row, l_hel,  &
                      HE_row_ptr, HE_val, l_he_val, HE_val_ptr, byrows )
      IF ( status /= 0 ) GO to 900
      byrows = .TRUE.
      WRITE( out, "( ' CALL CEH with byrows = .TRUE.' )" )
      CALL CEH( data, status, n, m, X, m, Y, HE_nel, HE_row, l_he_row, l_hel,  &
                      HE_row_ptr, HE_val, l_he_val, HE_val_ptr, byrows )
      IF ( status /= 0 ) GO to 900

!  compute the gradient and element Hessian values

      grlagf = .TRUE. ; byrows = .TRUE.
      WRITE( out, "( ' CALL CSGREH with grlagf = .TRUE. and byrows = .TRUE.')" )
      CALL CSGREH( data, status, n, m, X, grlagf, m, Y, J_ne, l_j, J_val,      &
                   J_col, J_row, HE_nel, HE_row, l_he_row, l_hel, HE_row_ptr,  &
                   HE_val, l_he_val, HE_val_ptr, byrows )
      IF ( status /= 0 ) GO to 900
      grlagf = .TRUE. ; byrows = .FALSE.
      WRITE( out, "(' CALL CSGREH with grlagf = .TRUE. and byrows = .FALSE.')" )
      CALL CSGREH( data, status, n, m, X, grlagf, m, Y, J_ne, l_j, J_val,      &
                   J_col, J_row, HE_nel, HE_row, l_he_row, l_hel, HE_row_ptr,  &
                   HE_val, l_he_val, HE_val_ptr, byrows )
      IF ( status /= 0 ) GO to 900
      grlagf = .FALSE. ; byrows = .TRUE.
      WRITE( out, "( ' CALL CSGREH with grlagf = .FALSE. and byrows = .TRUE.')")
      CALL CSGREH( data, status, n, m, X, grlagf, m, Y, J_ne, l_j, J_val,      &
                   J_col, J_row, HE_nel, HE_row, l_he_row, l_hel, HE_row_ptr,  &
                   HE_val, l_he_val, HE_val_ptr, byrows )
      IF ( status /= 0 ) GO to 900
      grlagf = .FALSE. ; byrows = .FALSE.
      WRITE( out, "(' CALL CSGREH with grlagf = .FALSE. and byrows = .FALSE.')")
      CALL CSGREH( data, status, n, m, X, grlagf, m, Y, J_ne, l_j, J_val,      &
                   J_col, J_row, HE_nel, HE_row, l_he_row, l_hel, HE_row_ptr,  &
                   HE_val, l_he_val, HE_val_ptr, byrows )
      IF ( status /= 0 ) GO to 900

!  compute a Hessian-vector product

      P = one
      goth = .FALSE.
      WRITE( out, "( ' Call CPROD with goth = .FALSE.' )" )
      CALL CPROD( data, status, n, m, goth, X, m, Y, P, HP )
      IF ( status /= 0 ) GO to 900
      goth = .TRUE.
      WRITE( out, "( ' Call CPROD with goth = .TRUE.' )" )
      CALL CPROD( data, status, n, m, goth, X, m, Y, P, HP )
      IF ( status /= 0 ) GO to 900

!  compute a Hessian-vector product ignoring the objective

      P = one
      goth = .FALSE.
      WRITE( out, "( ' Call CPROD1 with goth = .FALSE.' )" )
      CALL CPROD1( data, status, n, m, goth, X, m, Y, P, HP )
      IF ( status /= 0 ) GO to 900
      goth = .TRUE.
      WRITE( out, "( ' Call CPROD1 with goth = .TRUE.' )" )
      CALL CPROD1( data, status, n, m, goth, X, m, Y, P, HP )
      IF ( status /= 0 ) GO to 900

!  terminal exit

!     WRITE( out, "( ' CALL CREPRT' )" )
!     CALL CREPRT( CALLS, CPU )

      DEALLOCATE( X_type, H_row, H_col, HE_row, HE_row_ptr, HE_val_ptr, X,     &
                  X_l, X_u, G, Ji, Y, C_l, C_u, C, H_val, HE_val, H2_val,      &
                  H_band, J_val, J_row, J_col, J2_val, P, HP,                  &
                  X_names, C_names, EQUATION, LINEAR, stat = alloc_stat )
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

    END PROGRAM CUTEST_test_constrained_tools
