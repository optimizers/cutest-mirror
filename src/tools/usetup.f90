! THIS VERSION: CUTEST 1.0 - 23/12/2012 AT 15:10 GMT.

!-*-*-*-*-*-*-  C U T E S T    U S E T U P    S U B R O U T I N E  -*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 2003 version released in CUTEst, 22nd December 2012

      SUBROUTINE CUTEST_usetup( status, input, out, io_buffer, n, X, X_l, X_u )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      INTEGER, INTENT( IN ) :: input, out, io_buffer
      INTEGER, INTENT( INOUT ) :: n
      INTEGER, INTENT( OUT ) :: status
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: X, X_l, X_u

!  --------------------------------------------------------------
!  set up the input data for the unconstrained optimization tools
!  --------------------------------------------------------------

!  local variables

      INTEGER :: alloc_status
      CHARACTER ( LEN = 80 ) :: bad_alloc = REPEAT( ' ', 80 )

!  allocate space for the global workspace

      ALLOCATE( CUTEST_work_global( 1 ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'CUTEST_work_global' ; GO TO 910
      END IF

!  set the data

      CALL CUTEST_usetup_threadsafe( CUTEST_data_global,                       &
                                     CUTEST_work_global( 1 ),                  &
                                     status, input, out, io_buffer,            &
                                     n, X, X_l, X_u )
      CUTEST_data_global%threads = 1
      RETURN

!  allocation error

  910 CONTINUE
      status = 1
      IF ( out > 0 ) WRITE( out,                                               &
        "( /, ' ** SUBROUTINE USETUP: allocation error for ', A, ' status = ', &
       &  I0, /, ' Execution terminating ' )" ) TRIM( bad_alloc ), alloc_status
      RETURN

!  End of subroutine CUTEST_usetup

      END SUBROUTINE CUTEST_usetup

!-*-  C U T E S T    U S E T U P  _ t h r e a d e d   S U B R O U T I N E  -*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 2003 version released in CUTEst, 22nd December 2012

      SUBROUTINE CUTEST_usetup_threaded( status, input, out, threads,         &
                                         IO_BUFFERS, n, X, X_l, X_u )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      INTEGER, INTENT( IN ) :: input, out, threads
      INTEGER, INTENT( INOUT ) :: n
      INTEGER, INTENT( OUT ) :: status
      INTEGER, INTENT( IN ), DIMENSION( threads ) :: IO_BUFFERS
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: X, X_l, X_u

!  --------------------------------------------------------------
!  set up the input data for the unconstrained optimization tools
!  --------------------------------------------------------------

!  local variables

      INTEGER :: i, alloc_status
      CHARACTER ( LEN = 80 ) :: bad_alloc = REPEAT( ' ', 80 )

!  threads must be positive

      IF ( threads < 1 ) GO TO 940

!  allocate space for the threaded workspace data

      ALLOCATE( CUTEST_work_global( threads ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'CUTEST_work_global' ; GO TO 910
      END IF

!  set the data

      CALL CUTEST_usetup_threadsafe( CUTEST_data_global,                       &
                                     CUTEST_work_global( 1 ),                  &
                                     status, input, out, IO_BUFFERS( 1 ),      &
                                     n, X, X_l, X_u )
      CUTEST_data_global%threads = threads

!  initialize additional thread data

      DO i = 2, threads
        CALL CUTEST_initialize_thread( CUTEST_data_global,                     &
                                       CUTEST_work_global( i ), .FALSE.,       &
                                       status, alloc_status, bad_alloc )
        IF ( status /= 0 ) RETURN

!  ensure that each thread uses its own i/o buffer if required

        CUTEST_work_global( i )%io_buffer = IO_BUFFERS( i )

!  copy necessary starting addresses

        CUTEST_work_global( i )%ISTAJC( : CUTEST_data_global%n + 1 )           &
         = CUTEST_work_global( 1 )%ISTAJC( : CUTEST_data_global%n + 1 )
      END DO
      RETURN

!  allocation error

  910 CONTINUE
      status = 1
      IF ( out > 0 ) WRITE( out,                                               &
        "( /, ' ** SUBROUTINE USETUP: allocation error for ', A, ' status = ', &
       &  I0, /, ' Execution terminating ' )" ) TRIM( bad_alloc ), alloc_status
      RETURN

!  thread error

  940 CONTINUE
      status = 4
      IF ( out > 0 ) WRITE( out,                                               &
        "( /, ' ** SUBROUTINE USETUP: argument threads must be positive,',     &
       &  ' execution terminating ' )" )
      RETURN

!  End of subroutine CUTEST_usetup_threaded

      END SUBROUTINE CUTEST_usetup_threaded

!-*-  C U T E S T    U S E T U P _ t h r e a d s a f e  S U B R O U T I N E  -*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 77 version originally released in CUTE, 30th October, 1991
!   fortran 2003 version released in CUTEst, 4th November 2012

      SUBROUTINE CUTEST_usetup_threadsafe( data, work, status, input, out,     &
                                           io_buffer, n, X, X_l, X_u )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      TYPE ( CUTEST_data_type ), INTENT( INOUT ) :: data
      TYPE ( CUTEST_work_type ), INTENT( INOUT ) :: work
      INTEGER, INTENT( IN ) :: input, out, io_buffer
      INTEGER, INTENT( INOUT ) :: n
      INTEGER, INTENT( OUT ) :: status
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: X, X_l, X_u

!  --------------------------------------------------------------
!  set up the input data for the unconstrained optimization tools
!  --------------------------------------------------------------

!  local variables

      INTEGER :: i, ialgor, neltyp, ngrtyp, alloc_status
      LOGICAL :: debug
      INTEGER, PARAMETER :: lmin = 10000
      REAL ( KIND = wp ), DIMENSION( 2 ) :: OBFBND
      CHARACTER ( LEN = 8 ) :: pname
      CHARACTER ( LEN = 10 ) :: chtemp
      CHARACTER ( LEN = 80 ) :: bad_alloc = REPEAT( ' ', 80 )
      EXTERNAL :: RANGE

      CALL CPU_TIME( data%sutime )
      data%out = out
      work%io_buffer = io_buffer
      debug = .FALSE.
      debug = debug .AND. out > 0

!  input the problem dimensions

      READ( input, "( 10I8 )" ) data%n, data%ng, data%nel, data%ntotel,        &
        data%nvrels, data%nnza, data%ngpvlu, data%nepvlu, neltyp, ngrtyp
      n = data%n
      IF ( n <= 0 ) THEN
        CLOSE( input )
        IF ( out > 0 ) WRITE( out,                                             &
          "( /, ' ** SUNROUTINE USETUP: the problem uses no variables.',       &
         &       ' Execution terminating ' )" )
        status = 2 ; RETURN
      END IF
      IF ( data%ng <= 0 ) THEN
        CLOSE( input )
        IF ( out > 0 ) WRITE( out,                                             &
           "( /, ' ** SUBROUTINE USETUP: the problem is vacuous.',             &
          &      ' Execution terminating ' )" )
        status = 2 ; RETURN
      END IF
      IF ( SIZE( X ) < n ) THEN
        CLOSE( input )
        IF ( out > 0 ) WRITE( out, 2000 ) 'X', n
        status = 2 ; RETURN
      END IF
      IF ( SIZE( X_l ) < n ) THEN
        CLOSE( input )
        IF ( out > 0 ) WRITE( out, 2000 ) 'X_l', n
        status = 2 ; RETURN
      END IF
      IF ( SIZE( X_u ) < n ) THEN
        CLOSE( input )
        IF ( out > 0 ) WRITE( out, 2000 ) 'X_u', n
        status = 2 ; RETURN
      END IF

!  input the problem type

      READ( input, "( I2, A8 )" ) ialgor, pname

!  set useful integer values

      data%ng1 = data%ng + 1
      data%ngng = data%ng + data%ng
      data%nel1 = data%nel + 1

!  allocate integer workspace

      ALLOCATE( data%ISTADG( data%ng1 ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%ISTADG' ; GO TO 910
      END IF

      ALLOCATE( data%ISTGP( data%ng1 ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%ISTGP' ; GO TO 910
      END IF

      ALLOCATE( data%ISTADA( data%ng1 ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%ISTADA' ; GO TO 910
      END IF

      ALLOCATE( data%ISTAEV( data%nel1 ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%ISTAEV' ; GO TO 910
      END IF

      ALLOCATE( data%ISTEP( data%nel1 ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%ISTEP' ; GO TO 910
      END IF

      ALLOCATE( data%ITYPEG( data%ng), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%ITYPEG' ; GO TO 910
      END IF

      ALLOCATE( data%KNDOFC( data%ng ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%KNDOFC' ; GO TO 910
      END IF

      ALLOCATE( data%ITYPEE( data%nel ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%ITYPEE' ; GO TO 910
      END IF

      ALLOCATE( data%IELING( data%ntotel ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%IELING' ; GO TO 910
      END IF

      ALLOCATE( data%IELVAR( data%nvrels ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%IELVAR' ; GO TO 910
      END IF

      ALLOCATE( data%ICNA( data%nnza ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%ICNA' ; GO TO 910
      END IF

      ALLOCATE( data%ISTADH( data%nel1 ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%ISTADH' ; GO TO 910
      END IF

      ALLOCATE( data%INTVAR( data%nel1 ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%INTVAR' ; GO TO 910
      END IF

      ALLOCATE( data%IVAR( n ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%IVAR' ; GO TO 910
      END IF

      ALLOCATE( work%ICALCF( MAX( data%nel, data%ng ) ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'work%ICALCF' ; GO TO 910
      END IF

      ALLOCATE( data%ITYPEV( n ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%ITYPEV' ; GO TO 910
      END IF

      ALLOCATE( data%IWORK( 2 * n ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%IWORK' ; GO TO 910
      END IF

!  allocate real workspace

      ALLOCATE( data%A( data%nnza ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%A' ; GO TO 910
      END IF

      ALLOCATE( data%B( data%ng ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%B' ; GO TO 910
      END IF

      ALLOCATE( data%U( data%ng ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%U' ; GO TO 910
      END IF

      ALLOCATE( data%GPVALU( data%ngpvlu ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%GPVALU' ; GO TO 910
      END IF

      ALLOCATE( data%EPVALU( data%nepvlu ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%EPVALU' ; GO TO 910
      END IF

      ALLOCATE( data%ESCALE( data%ntotel ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%ESCALE' ; GO TO 910
      END IF

      ALLOCATE( data%GSCALE( data%ng ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%GSCALE' ; GO TO 910
      END IF

      ALLOCATE( work%GSCALE_used( data%ng ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'work%GSCALE_used' ; GO TO 910
      END IF

      ALLOCATE( data%VSCALE( n ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%VSCALE' ; GO TO 910
      END IF

      ALLOCATE( work%GVALS( data%ng, 3 ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'work%GVALS' ; GO TO 910
      END IF

      ALLOCATE( work%FT( data%ng ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'work%FT' ; GO TO 910
      END IF

!  allocate logical workspace

      ALLOCATE( data%INTREP( data%nel ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%INTREP' ; GO TO 910
      END IF

      ALLOCATE( data%GXEQX( data%ngng ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%GXEQX' ; GO TO 910
      END IF

!  allocate character workspace

      ALLOCATE( data%GNAMES( data%ng ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%GNAMES' ; GO TO 910
      END IF

      ALLOCATE( data%VNAMES( n ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%VNAMES' ; GO TO 910
      END IF

!  record the lengths of arrays

      data%ltypee = data%nel
      data%ltypeg = data%ng
      data%lstep  = data%nel1 
      data%lstgp = data%ng1
      data%lcalcf = MAX( data%nel, data%ng )
      data%lcalcg = data%ng
      data%lstaev = data%nel1
      data%lelvar = data%nvrels
      data%lntvar = data%nel1
      data%lstadh = data%nel1
      data%lvscal = n
      data%lepvlu = data%nepvlu
      data%lgpvlu = data%ngpvlu

!  assign initial guesses for allocatable array lengths

      work%lh_row = lmin
      work%lh_col = lmin
      work%lh_val = lmin

!  print out problem data. input the number of variables, groups, elements and 
!  the identity of the objective function group

      IF ( debug ) WRITE( out, 1100 ) pname, n, data%ng, data%nel
      data%pname = pname // '  '

!  input the starting addresses of the elements in each group, of the parameters
!  used for each group and of the nonzeros of the linear element in each group

      READ( input, 1010 ) ( data%ISTADG( i ), i = 1, data%ng1 )
      IF ( debug ) WRITE( out, 1110 ) 'ISTADG',                                &
        ( data%ISTADG( i ), i = 1, data%ng1 )
      READ( input, 1010 ) ( data%ISTGP( i ), i = 1, data%ng1 )
      IF ( debug ) WRITE( out, 1110 ) 'ISTGP ',                                &
        ( data%ISTGP( i ), i = 1, data%ng1 )
      READ( input, 1010 ) ( data%ISTADA( i ), i = 1, data%ng1 )
      IF ( debug ) WRITE( out, 1110 ) 'ISTADA',                                &
        ( data%ISTADA( i ), i = 1, data%ng1 )

!  input the starting addresses of the variables and parameters
!  in each element

      READ( input, 1010 ) ( data%ISTAEV( i ), i = 1, data%nel1 )
      IF ( debug ) WRITE( out, 1110 ) 'ISTAEV',                                &
        ( data%ISTAEV( i ), i = 1, data%nel1 )
      READ( input, 1010 ) ( data%ISTEP( i ), i = 1, data%nel1 )
      IF ( debug ) WRITE( out, 1110 ) 'ISTEP ',                                &
        ( data%ISTEP( i ), i = 1, data%nel1 )

!  input the group type of each group

      READ( input, 1010 ) ( data%ITYPEG( i ), i = 1, data%ng )
      IF ( debug ) WRITE( out, 1110 ) 'ITYPEG',                                &
        ( data%ITYPEG( i ), i = 1, data%ng )
      IF ( ialgor >= 2 ) THEN
        READ( input, 1010 )( data%KNDOFC( i ), i = 1, data%ng )
        IF ( debug ) WRITE( out, 1110 ) 'KNDOFC',                              &
          ( data%KNDOFC( i ), i = 1, data%ng )
        DO i = 1, data%ng
          IF ( ABS( data%KNDOFC( i ) ) >= 2 ) THEN
            CLOSE( input )
            IF ( out > 0 ) WRITE( out,                                         &
              "( /, ' ** Program USETUP: the problem includes general',        &
             &      ' constraints. Execution terminating ' )" )
            status = 2 ; RETURN
          END IF
        END DO
      END IF

!  input the element type of each element

      READ( input, 1010 ) ( data%ITYPEE( i ), i = 1, data%nel )
      IF ( debug ) WRITE( out, 1110 ) 'ITYPEE',                                &
        ( data%ITYPEE( i ), i = 1, data%nel )

!  input the number of internal variables for each element

      READ( input, 1010 ) ( data%INTVAR( i ), i = 1, data%nel )
      IF ( debug ) WRITE( out, 1110 ) 'INTVAR',                                &
        ( data%INTVAR( i ), i = 1, data%nel )

!  input the identity of each individual element

      READ( input, 1010 ) ( data%IELING( i ), i = 1, data%ntotel )
      IF ( debug ) WRITE( out, 1110 ) 'IELING',                                &
        ( data%IELING( i ), i = 1, data%ntotel )

!  input the variables in each group's elements

      data%nvrels = data%ISTAEV( data%nel1 ) - 1
      READ( input, 1010 ) ( data%IELVAR( i ), i = 1, data%nvrels )
      IF ( debug ) WRITE( out, 1110 ) 'IELVAR',                                &
        ( data%IELVAR( i ), i = 1, data%nvrels )

!  input the column addresses of the nonzeros in each linear element

      READ( input, 1010 ) ( data%ICNA( i ), i = 1, data%nnza )
      IF ( debug ) WRITE( out, 1110 ) 'ICNA  ',                                &
        ( data%ICNA( i ), i = 1, data%nnza )

!  input the values of the nonzeros in each linear element, the constant term 
!  in each group, the lower and upper bounds on the variables and the starting 
!  point for the minimization

      READ( input, 1020 ) ( data%A( i ), i = 1, data%nnza )
      IF ( debug ) WRITE( out, 1120 ) 'A     ',                                &
        ( data%A( i ), i = 1, data%nnza )
      READ( input, 1020 ) ( data%B( i ), i = 1, data%ng )
      IF ( debug ) WRITE( out, 1120 ) 'B     ',                                &
        ( data%B( i ), i = 1, data%ng )
      IF ( ialgor <= 2 ) THEN
        READ( input, 1020 ) ( X_l( i ), i = 1, n )
        IF ( debug ) WRITE( out, 1120 ) 'X_l    ', ( X_l( i ), i = 1, n )
        READ( input, 1020 ) ( X_u( i ), i = 1, n )
        IF ( debug ) WRITE( out, 1120 ) 'X_u    ', ( X_u( i ), i = 1, n )
      ELSE

!  use GVALS and FT as temporary storage for the constraint bounds

        READ( input, 1020 ) ( X_l( i ), i = 1, n ),                            &
          ( work%GVALS( i, 1 ), i = 1, data%ng )
        IF ( debug ) WRITE( out, 1120 ) 'X_l    ',                             &
          ( X_l( i ), i = 1, n ), ( work%GVALS( i, 1 ), i = 1, data%ng )
        READ( input, 1020 ) ( X_u( i ), i = 1, n ),                            &
          ( work%FT( i ), i = 1, data%ng )
        IF ( debug ) WRITE( out, 1120 ) 'X_u    ',                             &
          ( X_u( i ), i = 1, n ), ( work%FT( i ), i = 1, data%ng )
      END IF
      READ( input, 1020 ) ( X( i ), i = 1, n )
      IF ( debug ) WRITE( out, 1120 ) 'X     ', ( X( i ), i = 1, n )
      IF ( ialgor >= 2 ) THEN
        READ( input, 1020 )( data%U( i ), i = 1, data%ng )
        IF ( debug ) WRITE( out, 1120 ) 'U     ',                              &
          ( data%U( i ), i = 1, data%ng )
      END IF

!  input the parameters in each group

      READ( input, 1020 ) ( data%GPVALU( i ), i = 1, data%ngpvlu )
      IF ( debug ) WRITE( out, 1120 ) 'GPVALU',                                &
        ( data%GPVALU( i ), i = 1, data%ngpvlu )

!  input the parameters in each individual element

      READ( input, 1020 ) ( data%EPVALU( i ), i = 1, data%nepvlu )
      IF ( debug ) WRITE( out, 1120 ) 'EPVALU',                                &
        ( data%EPVALU( i ), i = 1, data%nepvlu )

!  input the scale factors for the nonlinear elements

      READ( input, 1020 ) ( data%ESCALE( i ), i = 1, data%ntotel )
      IF ( debug ) WRITE( out, 1120 ) 'ESCALE',                                &
        ( data%ESCALE( i ), i = 1, data%ntotel )

!  input the scale factors for the groups

      READ( input, 1020 ) ( data%GSCALE( i ), i = 1, data%ng )
      IF ( debug ) WRITE( out, 1120 ) 'GSCALE',                                &
        ( data%GSCALE( i ), i = 1, data%ng )

!  input the scale factors for the variables

      READ( input, 1020 ) ( data%VSCALE( i ), i = 1, n )
      IF ( debug ) WRITE( out, 1120 ) 'VSCALE',                                &
        ( data%VSCALE( i ), i = 1, n )

!  input the lower and upper bounds on the objective function

      READ( input, 1080 ) OBFBND( 1 ), OBFBND( 2 )
      IF ( debug ) WRITE( out, 1180 ) 'OBFBND', OBFBND( 1 ), OBFBND( 2 )

!  input a logical array which says whether an element has internal varaiables

      READ( input, 1030 ) ( data%INTREP( i ), i = 1, data%nel )
      IF ( debug ) WRITE( out, 1130 ) 'INTREP',                                &
        ( data%INTREP( i ), i = 1, data%nel )

!  input a logical array which says whether a group is trivial

      READ( input, 1030 ) ( data%GXEQX( i ), i = 1, data%ng )
      IF ( debug ) WRITE( out, 1130 ) 'GXEQX ',                                &
        ( data%GXEQX( i ), i = 1, data%ng )

!  input the names given to the groups and to the variables

      READ( input, 1040 ) ( data%GNAMES( i ), i = 1, data%ng )
      IF ( debug ) WRITE( out, 1140 ) 'GNAMES',                                &
        ( data%GNAMES( i ), i = 1, data%ng )
      READ( input, 1040 ) ( data%VNAMES( i ), i = 1, n )
      IF ( debug ) WRITE( out, 1140 ) 'VNAMES', ( data%VNAMES( i ), i = 1, n )

!  dummy input for the names given to the element and group types

      READ( input, 1040 ) ( chtemp, i = 1, neltyp )
      READ( input, 1040 ) ( chtemp, i = 1, ngrtyp )

!  input the type of each variable

      READ( input, 1010 ) ( data%ITYPEV( i ), i = 1, n )
      CLOSE( input )

!  partition the workspace arrays work%FUVALS, IWK and WK. Initialize certain 
!  portions of IWK

      work%firstg = .TRUE.
      data%ntotel = data%ISTADG( data%ng + 1 ) - 1
      data%nvrels = data%ISTAEV( data%nel + 1 ) - 1
      data%nnza = data%ISTADA( data%ng + 1 ) - 1

      CALL CUTEST_initialize_workspace(                                        &
             data%n, data%ng, data%nel, data%ntotel, data%nvrels, data%nnza,   &
             data%nvargp, data%IELING, data%ISTADG, data%IELVAR, data%ISTAEV,  &
             data%INTVAR, data%ISTADH, data%ICNA, data%ISTADA,                 &
             data%GXEQX, data%alllin, data%altriv, data%lfxi, data%lgxi,       &
             data%lhxi, data%lggfx, data%ldx, data%lgrjac, data%lnguvl,        &
             data%lnhuvl, data%ntotin, data%maxsel, data%maxsin, RANGE, 0,     &
             out, work%io_buffer, data%l_link_e_u_v, work%FUVALS, data%lfuval, &
             data%LINK_elem_uses_var, work%ISWKSP, work%ISTAJC, data%ISTAGV,   &
             data%ISVGRP, data%ISLGRP, data%IGCOLJ, data%ISYMMH, work%W_ws,    &
             work%W_el, work%W_in, work%H_el, work%H_in,                       &
             status, alloc_status, bad_alloc, work%array_status )
      IF ( status /= 0 ) RETURN

!  initialize the performance counters and variables

      work%nc2of = 0 ;  work%nc2og = 0 ; work%nc2oh = 0
      work%nc2cf = 0 ;  work%nc2cg = 0 ; work%nc2ch = 0 ; work%nhvpr = 0
      work%pnc = 0

      CALL CPU_TIME( data%sttime )
      data%sutime = data%sttime - data%sutime

      status = 0
      RETURN

  910 CONTINUE
      status = 1
      IF ( out > 0 ) WRITE( out,                                               &
        "( /, ' ** SUBROUTINE USETUP: allocation error for ', A, ' status = ', &
       &  I0, /, ' Execution terminating ' )" ) TRIM( bad_alloc ), alloc_status
      RETURN

!  non-executable statements

 1010 FORMAT( ( 10I8 ) )
 1020 FORMAT( ( 1P, 4D16.8 ) )
 1030 FORMAT( ( 72L1 ) )
 1040 FORMAT( ( 8A10 ) )
 1080 FORMAT( 1P, 2D16.8 )
 1100 FORMAT( A8, 3I8 )
 1110 FORMAT( 1X, A6, /, ( 1X, 10I8 ) )
 1120 FORMAT( 1X, A6, /, ( 1X, 1P, 4D16.8 ) )
 1130 FORMAT( 1X, A6, /, ( 1X, 72L1 ) )
 1140 FORMAT( 1X, A6, /, ( 1X, 8A10 ) )
 1180 FORMAT( 1X, A6, /, 1P, 2D16.6 )
 2000 FORMAT( /, ' ** SUBROUTINE USETUP: array length ', A, ' too small.', /, &
              ' -- Increase the dimension to at least ', I0, ' and restart.' )

!  End of subroutine CUTEST_usetup_threadsafe

      END SUBROUTINE CUTEST_usetup_threadsafe