! THIS VERSION: CUTEST 1.0 - 04/11/2012 AT 12:00 GMT.

!-*-*-*-*-*-*-  C U T E S T    C S E T U P    S U B R O U T I N E  -*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould with modifications by Ingrid Bongartz

!  History -
!   fortran 77 version originally released in CUTE, 30th October, 1991
!   fortran 2003 version released in CUTEst, 4th November 2012

      SUBROUTINE CSETUP( data, status, input, out, n, m, X, X_l, X_u, EQUATN, &
                         LINEAR, Y, C_l, C_u, efirst, lfirst, nvfrst )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      TYPE ( CUTEST_data_type ), INTENT( INOUT ) :: data
      INTEGER, INTENT( IN ) :: input, out
      INTEGER, INTENT( INOUT ) :: n, m
      INTEGER, INTENT( OUT ) :: status
      LOGICAL, INTENT( IN ) :: efirst, lfirst, nvfrst
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: X, X_l, X_u
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( m ) :: Y, C_l, C_u
      LOGICAL, INTENT( OUT ), DIMENSION( m ) :: EQUATN, LINEAR

!  ------------------------------------------------------------
!  set up the input data for the constrained optimization tools
!  ------------------------------------------------------------

!  local variables

      INTEGER :: ialgor, iprint, i, ig, j, jg, mend, meq, mlin, mmax
      INTEGER :: ii, k, iel, jwrk, kndv, nnlin, nend, neltyp, ngrtyp, itemp
      INTEGER :: alloc_status
      LOGICAL :: fdgrad, debug, ltemp
      CHARACTER ( LEN = 8 ) :: pname
      CHARACTER ( LEN = 10 ) :: ctemp
      CHARACTER ( LEN = 10 ) :: chtemp
      CHARACTER ( LEN = 80 ) :: bad_alloc = REPEAT( ' ', 80 )
      REAL ( KIND = wp ) :: atemp
      REAL ( KIND = wp ), PARAMETER :: zero = 0.0_wp
      REAL ( KIND = wp ), DIMENSION( 2 ) :: OBFBND
      EXTERNAL :: RANGE

      CALL CPU_TIME( data%sutime )
      data%out = out
      debug = .FALSE.
      debug = debug .AND. out > 0
      iprint = 0
      IF ( debug ) iprint = 3

!  input the problem dimensions

      READ( input, 1001 ) data%n, data%ng, data%nel, data%ntotel, data%nvrels, &
         data%nnza, data%ngpvlu, data%nepvlu, neltyp, ngrtyp
      n = data%n
      IF ( n <= 0 ) THEN
        CLOSE( input )
        IF ( out > 0 ) WRITE( out,                                             &
          "( /, ' ** SUBROUTINE CSETUP: the problem uses no variables.',       &
          &     ' Execution terminating ' )" )
        status = 2 ; RETURN
      END IF
      IF ( data%ng <= 0 ) THEN
        CLOSE( input )
        IF ( out > 0 ) WRITE( out,                                             &
          "( /, ' ** SUBROUTINE CSETUP: the problem is vacuous.',              &
          &     ' Execution terminating ' )" )
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

      READ( input, 1000 ) ialgor, pname

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

      ALLOCATE( data%ICALCF( MAX( data%nel, data%ng ) ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%ICALCF' ; GO TO 910
      END IF

      ALLOCATE( data%ITYPEV( n ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%ITYPEV' ; GO TO 910
      END IF

      ALLOCATE( data%IWORK( MAX( m, 2 * n ) ), STAT = alloc_status )
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

      ALLOCATE( data%VSCALE( n ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%VSCALE' ; GO TO 910
      END IF

      ALLOCATE( data%GVALS( data%ng, 3 ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%GVALS' ; GO TO 910
      END IF

      ALLOCATE( data%XT( n ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%XT' ; GO TO 910
      END IF

      ALLOCATE( data%DGRAD( n ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%DGRAD' ; GO TO 910
      END IF

      ALLOCATE( data%Q( n ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%Q' ; GO TO 910
      END IF

      ALLOCATE( data%FT( data%ng ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN
        bad_alloc = 'data%FT' ; GO TO 910
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

!     data%lstadg = MAX( 1, data%ng1 )
!     data%lstada = MAX( 1, data%ng1 )
!     data%lstaev = MAX( 1, data%nel1 )
!     data%lkndof = MAX( 1, data%ng )
!     data%leling = MAX( 1, data%ntotel )
!     data%lelvar = MAX( 1, data%nvrels )
!     data%licna = MAX( 1, data%nnza )
!     data%lstadh = MAX( 1, data%nel1 )
!     data%lntvar = MAX( 1, data%nel1 )
!     data%lcalcf = MAX( 1, data%nel, data%ng )
!     data%lcalcg = MAX( 1, data%ng )
!     data%la = MAX( 1, data%nnza )
!     data%lb = MAX( 1, data%ng )
!     data%lu = MAX( 1, data%ng )
!     data%lescal = MAX( 1, data%ntotel )
!     data%lgscal = MAX( 1, data%ng )
!     data%lvscal = MAX( 1, n )
!     data%lft = MAX( 1, data%ng )
!     data%lgvals = MAX( 1, data%ng )
!     data%lintre = MAX( 1, data%nel )
!     data%lgxeqx = MAX( 1, data%ngng )
!     data%lgpvlu = MAX( 1, data%ngpvlu )
!     data%lepvlu = MAX( 1, data%nepvlu )
!     LSTGP = MAX( 1, data%ng1 )
!     LSTEP = MAX( 1, data%nel1 )
!     LTYPEG = MAX( 1, data%ng )
!     LTYPEE = MAX( 1, data%nel )
!     LIVAR = MAX( 1, n )
!     LBL = MAX( 1, n )
!     LBU = MAX( 1, n )
!     LX = MAX( 1, n )
!     LXT = MAX( 1, n )
!     LDGRAD = MAX( 1, n )
!     LQ = MAX( 1, n )

!  print out problem data. input the number of variables, groups, elements and 
!  the identity of the objective function group

!     IF ( ialgor == 2 ) THEN
!       READ( input, 1002 ) nslack, data%nobjgr
!     ELSE
!       nslack = 0
!     END IF
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

!  input the starting addresses of the variables and parameters in each element

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
         READ( input, 1010 ) ( data%KNDOFC( i ), i = 1, data%ng )
         IF ( debug ) WRITE( out, 1110 ) 'KNDOFC',                             &
        ( data%KNDOFC( i ), i = 1, data%ng )
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
!  in each group, the lower and upper bounds on the variables.

      READ( input, 1020 ) ( data%A( i ), i = 1, data%nnza )
      IF ( debug ) WRITE( out, 1120 ) 'A     ', ( data%A( i ), i = 1, data%nnza)
      READ( input, 1020 ) ( data%B( i ), i = 1, data%ng )
      IF ( debug ) WRITE( out, 1120 ) 'B     ', ( data%B( i ), i = 1, data%ng )
      IF ( ialgor <= 2 ) THEN
         READ( input, 1020 ) ( X_l( i ), i = 1, n )
         IF ( debug ) WRITE( out, 1120 ) 'X_l    ', ( X_l( i ), i = 1, n )
         READ( input, 1020 ) ( X_u( i ), i = 1, n )
         IF ( debug ) WRITE( out, 1120 ) 'X_u    ', ( X_u( i ), i = 1, n )
      ELSE

!  use GVALS and FT as temporary storage for the constraint bounds

         READ( input, 1020 ) ( X_l( i ), i = 1, n ),                           &
           ( data%GVALS( i, 1 ), i = 1, data%ng )
         IF ( debug ) WRITE( out, 1120 ) 'X_l    ',                            &
           ( X_l( i ), i = 1, n ), ( data%GVALS( i, 1 ), i = 1, data%ng )
         READ( input, 1020 ) ( X_u( i ), i = 1, n ),                           &
           ( data%FT( i ), i = 1, data%ng )
         IF ( debug ) WRITE( out, 1120 ) 'X_u    ',                            &
           ( X_u( i ), i = 1, n ), ( data%FT( i ), i = 1, data%ng )
      END IF

!   input the starting point for the minimization

      READ( input, 1020 ) ( X( i ), i = 1, n )
      IF ( debug ) WRITE( out, 1120 ) 'X     ', ( X( i ), i = 1, n )
      IF ( ialgor >= 2 ) THEN
         READ( input, 1020 )( data%U( i ), i = 1, data%ng )
         IF ( debug ) WRITE( out, 1120 ) 'U     ',                             &
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

!  input a logical array which says whether an element has internal variables

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
      IF ( debug ) WRITE( out, 1140 ) 'VNAMES',                                &
        ( data%VNAMES( i ), i = 1, n )

!  dummy input for the names given to the element and group types

      READ( input, 1040 ) ( CHTEMP, i = 1, neltyp )
      READ( input, 1040 ) ( CHTEMP, i = 1, ngrtyp )

!  input the type of each variable.

      READ( input, 1010 ) ( data%ITYPEV( i ), i = 1, n )

!  consider which groups are constraints. Of these, decide which are
!  equations, which are linear, allocate starting values for the
!  Lagrange multipliers and set lower and upper bounds on any inequality 
!  constraints. Reset KNDOFC to point to the list of constraint groups.

      mmax = MIN( SIZE( LINEAR ), SIZE( EQUATN ),                              &
                  SIZE( C_l ), SIZE( C_u ), SIZE( Y )  )
      m = 0
      DO i = 1, data%ng
        IF ( data%KNDOFC( i ) == 1 ) THEN
          data%KNDOFC( i ) = 0
        ELSE
          m = m + 1
          IF ( m <= mmax ) THEN
            Y( m ) = data%U( i )
            LINEAR( m ) =                                                      &
              data%GXEQX( i ) .AND. data%ISTADG( i ) >= data%ISTADG( i + 1 )
            IF ( data%KNDOFC( i ) == 2 ) THEN
              EQUATN( m ) = .TRUE.
              C_l( m ) = zero ; C_u( m ) = zero
            ELSE
              EQUATN( m ) = .FALSE.
              C_l( m ) = data%GVALS( i, 1 ) ; C_u( m ) = data%FT( i )
            END IF
          END IF
          data%KNDOFC( i ) = m
        END IF
      END DO
      IF ( m == 0 .AND. out > 0 ) WRITE( out,                                  &
        "( /, ' ** SUBROUTINE CSETUP: ** Warning. The problem has',            &
       &      ' no general constraints. ', /,                                  &
       &      ' Other tools may be preferable' )" )
      IF ( SIZE( Y ) < m ) THEN
        CLOSE( input )
        IF ( out > 0 ) WRITE( out, 2000 ) 'V', m
        status = 2 ; RETURN
      END IF
      IF ( SIZE( C_l ) < m ) THEN
        CLOSE( input )
        IF ( out > 0 ) WRITE( out, 2000 ) 'C_l', m
        status = 2 ; RETURN
      END IF
      IF ( SIZE( C_u ) < m ) THEN
        CLOSE( input )
        IF ( out > 0 ) WRITE( out, 2000 ) 'C_u', m
        status = 2 ; RETURN
      END IF
      IF ( SIZE( EQUATN ) < m ) THEN
        CLOSE( input )
        IF ( out > 0 ) WRITE( out, 2000 ) 'EQUATN', m
        status = 2 ; RETURN
      END IF
      IF ( SIZE( LINEAR ) < m ) THEN
        CLOSE( input )
        IF ( out > 0 ) WRITE( out, 2000 ) 'LINEAR', m
        status = 2 ; RETURN
      END IF
      data%numvar = n
      data%numcon = m

      IF ( nvfrst ) THEN

!  set addresses in IWORK to reorder variables

        kndv = 0 
        jwrk = kndv + n

!  initialize jwrk and kndv

        DO j = 1, n
          data%IWORK( kndv + j ) = 0
          data%IWORK( jwrk + j ) = j
        END DO

!  now identify and count nonlinear variables; keep separate counts for 
!  nonlinear objective and Jacobian variables.
!  data%IWORK(kndv + j) = 0 ==> j linear everywhere
!  data%IWORK(kndv + j) = 1 ==> j linear in objective, nonlinear in constraints
!  data%IWORK(kndv + j) = 2 ==> j linear in constraints, nonlinear in objective
!  data%IWORK(kndv + j) = 3 ==> j nonlinear everywhere

        nnlin = 0
        data%nnov = 0
        data%nnjv = 0
        DO ig = 1, data%ng
          i = data%KNDOFC( ig )
          DO  ii = data%ISTADG( ig ), data%ISTADG( ig + 1 ) - 1
            iel = data%IELING( ii )
            DO k = data%ISTAEV( iel ), data%ISTAEV( iel + 1 ) - 1
              j = data%IELVAR( k )
              IF ( i > 0 ) THEN
                IF ( data%IWORK( kndv + j ) == 0 ) THEN
                  data%IWORK( kndv + j ) = 1
                  data%nnjv = data%nnjv + 1
                  nnlin = nnlin + 1
                ELSE IF ( data%IWORK( kndv + j ) == 2 ) THEN
                  data%IWORK( kndv + j ) = 3
                  data%nnjv = data%nnjv + 1
                END IF
              ELSE
                IF ( data%IWORK( kndv + j ) == 0 ) THEN
                  data%IWORK( kndv + j ) = 2
                  data%nnov = data%nnov + 1
                  nnlin = nnlin + 1
                ELSE IF ( data%IWORK( kndv + j ) == 1 ) THEN
                  data%IWORK( kndv + j ) = 3
                  data%nnov = data%nnov + 1
                END IF
              END IF
            END DO
          END DO
          IF ( .NOT. data%GXEQX( ig ) ) THEN
            DO ii = data%ISTADA( ig ), data%ISTADA( ig + 1 ) - 1
              j = data%ICNA( ii )
              IF ( i > 0 ) THEN
                IF ( data%IWORK( kndv + j ) == 0 ) THEN
                  data%IWORK( kndv + j ) = 1
                  data%nnjv = data%nnjv + 1
                  nnlin = nnlin + 1
                ELSE IF ( data%IWORK( kndv + j ) == 2 ) THEN
                  data%IWORK( kndv + j ) = 3
                  data%nnjv = data%nnjv + 1
                END IF
              ELSE
                IF ( data%IWORK( kndv + j ) == 0 ) THEN
                  data%IWORK( kndv + j ) = 2
                  data%nnov = data%nnov + 1
                  nnlin = nnlin + 1
                ELSE IF ( data%IWORK( kndv + j ) == 1 ) THEN
                  data%IWORK( kndv + j ) = 3
                  data%nnov = data%nnov + 1
                END IF
              END IF
            END DO
          END IF
        END DO
        IF ( nnlin == 0 .OR. ( data%nnov == n .AND. data%nnjv == n ) ) GO TO 600
        IF ( nnlin == n ) GO TO 500

!  reorder the variables so that all nonlinear variables occur before the
!  linear ones

        nend = n

!  run forward through the variables until a linear variable is encountered

        DO 420 i = 1, n
          IF ( i > nend ) GO TO 430
          IF ( data%IWORK( kndv + i ) == 0 ) THEN

!  variable i is linear. Now, run backwards through the variables until a 
!  nonlinear one is encountered

            DO j = nend, i, - 1
              IF ( data%IWORK( kndv + j ) > 0 ) THEN 
                nend = j - 1

!  interchange the data for variables i and j

                itemp = data%IWORK( jwrk + i )
                data%IWORK( jwrk + i ) = data%IWORK( jwrk + j )
                data%IWORK( jwrk + j ) = itemp
                itemp = data%IWORK( kndv + i )
                data%IWORK( kndv + i ) = data%IWORK( kndv + j )
                data%IWORK( kndv + j ) = itemp
                atemp = X_l( i ) ; X_l( i ) = X_l( j ) ; X_l( j ) = atemp
                atemp = X_u( i ) ; X_u( i ) = X_u( j ) ; X_u( j ) = atemp
                atemp = X( i ) ; X( i ) = X( j ) ; X( j ) = atemp
                atemp = data%VSCALE( i )
                data%VSCALE( i ) = data%VSCALE( j )
                data%VSCALE( j ) = atemp
                ctemp = data%VNAMES( i )
                data%VNAMES( i ) = data%VNAMES( j )
                data%VNAMES( j ) = ctemp
                GO TO 420
              END IF
            END DO
            GO TO 430
          END IF
  420   CONTINUE
  430   CONTINUE 

!  change entries in IELVAR and ICNA to reflect reordering of variables

        DO i = 1, data%nvrels
          j = data%IELVAR( i )
          data%IELVAR( i ) = data%IWORK( jwrk + j ) 
        END DO
        DO i = 1, data%nnza
          j = data%ICNA( i )
          data%ICNA( i ) = data%IWORK( jwrk + j )
        END DO
        DO j = 1, n
           data%IWORK( jwrk + j ) = j
        END DO
  500   CONTINUE
        IF ( ( data%nnov == nnlin .AND. data%nnjv == nnlin )                   &
           .OR. ( data%nnov == 0 ) .OR. ( data%nnjv == 0 ) ) GO TO 600

!  reorder the nonlinear variables so that the smaller set (nonlinear objective 
!  or nonlinear Jacobian) occurs at the beginning of the larger set

        nend = nnlin
        IF ( data%nnjv <= data%nnov ) THEN

!  put the nonlinear Jacobian variables first. Reset data%nnov to indicate all 
!  nonlinear variables are treated as  nonlinear objective variables.

          data%nnov = nnlin
          DO 520 i = 1, nnlin 
            IF ( i > nend ) GO TO 530
            IF ( data%IWORK( kndv + i ) == 2 ) THEN

!  variable i is linear in the Jacobian. Now, run backwards through the 
!  variables until a nonlinear Jacobian variable is encountered

              DO j = nend, i, - 1
                IF ( data%IWORK( kndv + j ) == 1 .OR.                          &
                     data%IWORK( kndv + j ) == 3 ) THEN 
                  nend = j - 1

!  Interchange the data for variables i and j

                  itemp = data%IWORK( jwrk + i )
                  data%IWORK( jwrk + i ) = data%IWORK( jwrk + j )
                  data%IWORK( jwrk + j ) = itemp
                  itemp = data%IWORK( kndv + i )
                  data%IWORK( kndv + i ) = data%IWORK( kndv + j )
                  data%IWORK( kndv + j ) = itemp
                  atemp = X_l( i ) ; X_l( i ) = X_l( j ) ; X_l( j ) = atemp
                  atemp = X_u( i ) ; X_u( i ) = X_u( j ) ; X_u( j ) = atemp
                  atemp = X( i ) ; X( i ) = X( j ) ; X( j ) = atemp
                  atemp = data%VSCALE( i )
                  data%VSCALE( i ) = data%VSCALE( j )
                  data%VSCALE( j ) = atemp
                  ctemp = data%VNAMES( i )
                  data%VNAMES( i ) = data%VNAMES( j )
                  data%VNAMES( j ) = ctemp
                  GO TO 520
                END IF
              END DO
              GO TO 530
            END IF
  520     CONTINUE
  530     CONTINUE
        ELSE

!  put the nonlinear objective variables first. Reset data%nnjv to indicate all
!  nonlinear variables are treated as nonlinear Jacobian variables.

          data%nnjv = nnlin
          DO 550 i = 1, nnlin 
            IF ( i > nend ) GO TO 560
            IF ( data%IWORK( kndv + i ) == 1 ) THEN

!  variable i is linear in the objective. Now, run backwards through the 
!  variables until a nonlinear objective variable is encountered

              DO 540 j = nend, i, - 1
                 IF ( data%IWORK( kndv + j ) > 1 ) THEN 
                   nend = j - 1

!  interchange the data for variables i and j

                   itemp = data%IWORK( jwrk + i )
                   data%IWORK( jwrk + i ) = data%IWORK( jwrk + j )
                   data%IWORK( jwrk + j ) = itemp
                   itemp = data%IWORK( kndv + i )
                   data%IWORK( kndv + i ) = data%IWORK( kndv + j )
                   data%IWORK( kndv + j ) = itemp
                   atemp = X_l( i ) ;  X_l( i ) = X_l( j ) ;  X_l( j ) = atemp
                   atemp = X_u( i ) ;  X_u( i ) = X_u( j ) ;  X_u( j ) = atemp
                   atemp = X( i ) ;  X( i ) = X( j ) ;  X( j ) = atemp
                   atemp = data%VSCALE( i )
                   data%VSCALE( i ) = data%VSCALE( j )
                   data%VSCALE( j ) = atemp
                   ctemp = data%VNAMES( i )
                   data%VNAMES( i ) = data%VNAMES( j )
                   data%VNAMES( j ) = ctemp
                   GO TO 550
                 END IF
  540         CONTINUE
              GO TO 560
            END IF
  550     CONTINUE
  560     CONTINUE
        END IF

!  change entries in IELVAR and ICNA to reflect reordering of variables

        DO i = 1, data%nvrels
          j = data%IELVAR( i )
          data%IELVAR( i ) = data%IWORK( jwrk + j ) 
        END DO
        DO i = 1, data%nnza
          j = data%ICNA( i )
          data%ICNA( i ) = data%IWORK( jwrk + j )
        END DO
  600   CONTINUE
      END IF

!  allocate and initialize workspace

      data%firstg = .TRUE.
      fdgrad = .FALSE.

      data%ntotel = data%ISTADG( data%ng + 1 ) - 1
      data%nvrels = data%ISTAEV( data%nel + 1 ) - 1
      data%nnza = data%ISTADA( data%ng + 1 ) - 1

      CALL CUTEST_initialize_workspace(                                        &
             data%n, data%ng, data%nel,                                        &
             data%ntotel, data%nvrels, data%nnza, data%n,                      &
             data%nvargp, data%IELING, data%ISTADG, data%IELVAR, data%ISTAEV,  &
             data%INTVAR, data%ISTADH, data%ICNA, data%ISTADA, data%ITYPEE,    &
             data%GXEQX, data%INTREP, data%alllin, data%altriv, .FALSE.,       &
             .FALSE., data%lfxi, data%lgxi, data%lhxi,                         &
             data%lggfx, data%ldx, data%lnguvl, data%lnhuvl,                   &
             data%ntotin, data%ntype, data%nsets, data%maxsel,                 &
             RANGE, 0, out, data%io_buffer,                                    &
!  workspace
             data%lwtran, data%litran, data%lwtran_min, data%litran_min,       &
             data%l_link_e_u_v, data%llink_min, data%FUVALS, data%lfuval,      &
             data%ITRANS, data%LINK_elem_uses_var, data%WTRANS,                &
             data%ISYMMD, data%ISWKSP, data%ISTAJC, data%ISTAGV,               &
             data%ISVGRP, data%ISLGRP, data%IGCOLJ, data%IVALJR,               &
             data%IUSED, data%ITYPER, data%ISSWTR, data%ISSITR,                &
             data%ISET, data%ISVSET, data%INVSET, data%LIST_elements,          &
             data%ISYMMH, data%IW_asmbl, data%NZ_comp_w, data%W_ws,            &
             data%W_el, data%W_in, data%H_el, data%H_in,                       &
             status, alloc_status, bad_alloc, data%skipg, KNDOFG = data%KNDOFC )

!     CALL INITW( n, data%ng, data%nel, data%IELING, data%leling,              &
!         data%ISTADG, data%lstadg, data%IELVAR, data%lelvar, data%ISTAEV,     &
!         data%lstaev, data%INTVAR, data%lntvar, data%ISTADH, data%lstadh,     &
!         data%ICNA, data%licna, data%ISTADA, data%lstada, data%ITYPEE,        &
!         data%lintre, data%GXEQX, data%lgxeqx, data%INTREP, data%lintre,      &
!         lfuval, data%altriv, .TRUE., FDGRAD, data%lfxi, LGXI, LHXI, LGGFX,   &
!         data%ldx, data%lgrjac, data%lqgrad, data%lbreak, data%lp, data%lxcp, &
!         data%lx0, data%lgx0, data%ldeltx, data%lbnd, data%lwkstr,            &
!         data%lsptrs, data%lselts, data%lindex, data%lswksp, data%lstagv,     &
!         data%lstajc, data%liused, data%lfreec, data%lnnonz, data%lnonz2,     &
!         data%lsymmd, data%lsymmh, data%lslgrp, data%lsvgrp, data%lgcolj,     &
!         data%lvaljr, data%lsend,  data%lnptrs, data%lnelts, data%lnndex,     &
!         data%lnwksp, data%lnstgv, data%lnstjc, data%lniuse,  data%lnfrec,    &
!         data%lnnnon, data%lnnno2, data%lnsymd, data%lnsymh, data%lnlgrp,     &
!         data%lnvgrp, data%lngclj, data%lnvljr, data%lnqgrd, data%lnbrak,     &
!         data%lnp, data%lnbnd, data%lnfxi,  data%lngxi,  data%lnguvl,         &
!         data%lnhxi,  data%lnhuvl, data%lnggfx, data%lndx, data%lngrjc,       &
!         data%liwk2, data%lwk2, data%maxsin, data%ninvar, data%ntype,         &
!         data%nsets, data%maxsel, data%lstype, data%lsswtr, data%lssiwt,      &
!         data%lsiwtr, data%lswtra, data%lntype, data%lnswtr, data%lnsiwt,     &
!         data%lniwtr, data%lnwtra, data%lsiset, data%lssvse, data%lniset,     &
!         data%lnsvse, RANGE, data%IWORK(IWRK + 1), liwork, data%WRK, lwork,   &
!         iprint, out, status )
      IF ( status /= 0 ) RETURN

!  shift the starting addresses for the real workspace relative to WRK

!     data%lqgrad = data%lqgrad + WRK
!     data%lbreak = data%lbreak + WRK
!     data%lp = data%lp + WRK
!     data%lxcp = data%lxcp + WRK
!     data%lx0 = data%lx0 + WRK
!     data%lgx0 = data%lgx0 + WRK
!     data%ldeltx = data%ldeltx + WRK
!     data%lbnd = data%lbnd + WRK
!     data%lswtra = data%lswtra + WRK
!     data%lwkstr = data%lwkstr + WRK

!  shift the starting addresses for the integer workspace relative to IWRK

!     data%lsptrs = data%lsptrs + IWRK
!     data%lselts = data%lselts + IWRK
!     data%lindex = data%lindex + IWRK
!     data%lswksp = data%lswksp + IWRK
!     data%lstagv = data%lstagv + IWRK
!     data%lstajc = data%lstajc + IWRK
!     data%liused = data%liused + IWRK
!     data%lfreec = data%lfreec + IWRK
!     data%lnnonz = data%lnnonz + IWRK
!     data%lnonz2 = data%lnonz2 + IWRK
!     data%lsymmd = data%lsymmd + IWRK
!     data%lsymmh = data%lsymmh + IWRK
!     data%lslgrp = data%lslgrp + IWRK
!     data%lsvgrp = data%lsvgrp + IWRK
!     data%lgcolj = data%lgcolj + IWRK
!     data%lvaljr = data%lvaljr + IWRK
!     data%lstype = data%lstype + IWRK
!     data%lsswtr = data%lsswtr + IWRK
!     data%lssiwt = data%lssiwt + IWRK
!     data%lsiwtr = data%lsiwtr + IWRK
!     data%lsiset = data%lsiset + IWRK
!     data%lssvse = data%lssvse + IWRK
!     data%lsend = data%lsend + IWRK
      IF ( .NOT. ( efirst .OR. lfirst ) .OR. m == 0 ) GO TO 340

!  record which group is associated with each constraint

      meq = 0
      mlin = 0
      DO ig = 1, data%ng
        i = data%KNDOFC( ig )
        IF ( i > 0 ) THEN
          data%IWORK( i ) = ig
          IF ( EQUATN( i ) ) meq = meq + 1
          IF ( LINEAR( i ) ) mlin = mlin + 1
        END IF
      END DO
      IF ( lfirst ) THEN
        IF ( mlin == 0 .OR. mlin == m ) GO TO 130

!  reorder the constraints so that the linear constraints occur before the
!  nonlinear ones

        mend = m

!  run forward through the constraints until a nonlinear constraint is found

        DO 120 i = 1, m
          IF ( i > mend ) GO TO 130
          ig = data%IWORK( i )
!         WRITE(6,*) ' group ', ig, ' type ', i, ' equal? ', EQUATN( i )
          IF ( .NOT. LINEAR( i ) ) THEN

!  constraint i is nonlinear. Now, run backwards through the constraints until 
!  a linear one is encountered

            DO j = mend, i, - 1
              jg = data%IWORK( j )
!             write(6,*) ' group ', jg, ' type ', j, ' linear? ', LINEAR( j )
              IF ( LINEAR( j ) ) THEN
!               write(6,*) ' swaping constraints ', i, ' and ', j
                 mend = j - 1

!  interchange the data for constraints i and j

                 data%IWORK( i ) = jg
                 data%IWORK( j ) = ig
                 data%KNDOFC( ig ) = j
                 data%KNDOFC( jg ) = i
                 ltemp = LINEAR( i )
                 LINEAR( i ) = LINEAR( j )
                 LINEAR( j ) = ltemp
                 ltemp = EQUATN( i )
                 EQUATN( i ) = EQUATN( j )
                 EQUATN( j ) = ltemp
                 atemp = Y( i ) ; Y( i ) = Y( j ) ; Y( j ) = atemp
                 atemp = C_l( i ) ; C_l( i ) = C_l( j ) ; C_l( j ) = atemp
                 atemp = C_u( i ) ; C_u( i ) = C_u( j ) ; C_u( j ) = atemp
                 GO TO 120
               END IF
            END DO
            GO TO 130
          END IF
  120   CONTINUE
  130   CONTINUE
        IF ( efirst ) THEN
          IF ( meq == 0 .OR. meq == m ) GO TO 260

!  reorder the linear constraints so that the equations occur before the 
!  inequalities

          mend = mlin
          DO 220 i = 1, mlin
            IF ( i > mend ) GO TO 230
            ig = data%IWORK( i )
!              write(6,*) ' group ', ig, ' type ', i, ' equation? ',           &
!                           EQUATN( i )
            IF ( .NOT. EQUATN( i ) ) THEN

!  constraint i is an inequality. Now, run backwards through the constraints 
!  until an equation is encountered

              DO j = mend, i, - 1
                jg = data%IWORK( j )
!               write(6,*) ' group ', jg, ' type ', j, ' equation? ', EQUATN( j)
                IF ( EQUATN( j ) ) THEN
!                 write(6,*) ' swaping constraints ', i,                       &
!                            ' and ', j
                  mend = j - 1

!  interchange the data for constraints i and j

                  data%IWORK( i ) = jg
                  data%IWORK( j ) = ig
                  data%KNDOFC( ig ) = j
                  data%KNDOFC( jg ) = i
                  ltemp = LINEAR( i )
                  LINEAR( i ) = LINEAR( j )
                  LINEAR( j ) = ltemp
                  ltemp = EQUATN( i )
                  EQUATN( i ) = EQUATN( j )
                  EQUATN( j ) = ltemp
                  atemp = Y( i ) ; Y( i ) = Y( j ) ; Y( j ) = atemp
                  atemp = C_l( i ) ; C_l( i ) = C_l( j ) ; C_l( j ) = atemp
                  atemp = C_u( i ) ; C_u( i ) = C_u( j ) ; C_u( j ) = atemp
                  GO TO 220
                END IF
              END DO
              GO TO 230
            END IF
  220     CONTINUE
  230     CONTINUE

!  reorder the nonlinear constraints so that the equations occur  before the 
!  inequalities

          mend = m
          DO 250 i = mlin + 1, m
            IF ( i > mend ) GO TO 260
            ig = data%IWORK( i )
!           WRITE( 6, * ) ' group ', ig, ' type ', i, ' equation? ', EQUATN( i )
            IF ( .NOT. EQUATN( i ) ) THEN

!  constraint i is an inequality. Now, run backwards through the constraints 
!  until an equation is encountered

              DO j = mend, i, - 1
                jg = data%IWORK( j )
!               write(6,*) ' group ', jg, ' type ', j, ' equation? ', EQUATN( j)
                IF ( EQUATN( j ) ) THEN
!                 write(6,*) ' swaping constraints ', i, ' and ', j
                  mend = j - 1

!  interchange the data for constraints i and j

                  data%IWORK( i ) = jg
                  data%IWORK( j ) = ig
                  data%KNDOFC( ig ) = j
                  data%KNDOFC( jg ) = i
                  ltemp = LINEAR( i )
                  LINEAR( i ) = LINEAR( j )
                  LINEAR( j ) = ltemp
                  ltemp = EQUATN( i )
                  EQUATN( i ) = EQUATN( j )
                  EQUATN( j ) = ltemp
                  atemp = Y( i ) ; Y( i ) = Y( j ) ; Y( j ) = atemp
                  atemp = C_l( i ) ; C_l( i ) = C_l( j ) ; C_l( j ) = atemp
                  atemp = C_u( i ) ; C_u( i ) = C_u( j ) ; C_u( j ) = atemp
                  GO TO 250
                END IF
              END DO
              GO TO 260
            END IF
  250     CONTINUE
  260     CONTINUE
        END IF
      ELSE
        IF ( efirst ) THEN
          IF ( meq == 0 .OR. meq == m ) GO TO 330

!  reorder the constraints so that the equations occur before the inequalities

          mend = m
          DO 320 i = 1, m
            IF ( i > mend ) GO TO 330
            ig = data%IWORK( i )
!              WRITE(6,*) ' group ', ig, ' type ', i, ' equation? ', EQUATN( i )
            IF ( .NOT. EQUATN( i ) ) THEN

!  constraint i is an inequality. Now, run backwards through the constraints 
!  until an equation is encountered

              DO j = mend, i, - 1
                jg = data%IWORK( j )
!               write(6,*) ' group ', jg, ' type ', j,                         &
!                          ' equation? ', EQUATN( j )
                IF ( EQUATN( j ) ) THEN
!                 write(6,*) ' swaping constraints ', i,' and ', j
                  mend = j - 1

!  interchange the data for constraints i and j

                  data%IWORK( i ) = jg
                  data%IWORK( j ) = ig
                  data%KNDOFC( ig ) = j
                  data%KNDOFC( jg ) = i
                  ltemp = LINEAR( i )
                  LINEAR( i ) = LINEAR( j )
                  LINEAR( j ) = ltemp
                  ltemp = EQUATN( i )
                  EQUATN( i ) = EQUATN( j )
                  EQUATN( j ) = ltemp
                  atemp = Y( i ) ; Y( i ) = Y( j ) ; Y( j ) = atemp
                  atemp = C_l( i ) ; C_l( i ) = C_l( j ) ; C_l( j ) = atemp
                  atemp = C_u( i ) ; C_u( i ) = C_u( j ) ; C_u( j ) = atemp
                  GO TO 320
                END IF
              END DO
              GO TO 330
            END IF
  320     CONTINUE
  330     CONTINUE
        END IF
      END IF

!  initialize the performance counters and variables

 340  CONTINUE
      data%nc2of = 0
      data%nc2og = 0
      data%nc2oh = 0
      data%nc2cf = 0
      data%nc2cg = 0
      data%nc2ch = 0
      data%nhvpr = 0
      data%pnc = m

      CALL CPU_TIME( data%sttime )
      data%sutime = data%sttime - data%sutime

      status = 0
      RETURN

  910 CONTINUE
      status = 1
      IF ( out > 0 ) WRITE( out,                                               &
        "( /, ' ** SUBROUTINE CSETUP: allocation error for ', A, ' status = ', &
       &  I0, /, ' Execution terminating ' )" ) bad_alloc, alloc_status
      RETURN

!  non-executable statements

 1000 FORMAT( I2, A8 )
 1001 FORMAT( 10I8 )
!1002 FORMAT( 2I8 )
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
 2000 FORMAT( /, ' ** SUBROUTINE CSETUP: array length ', A, ' too small.', /, &
              ' -- Increase the dimension to at least ', I0, ' and restart.' )

!  End of subroutine CSETUP

      END SUBROUTINE CSETUP
