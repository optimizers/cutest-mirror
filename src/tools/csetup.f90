      SUBROUTINE CSETUP( data, INPUT, IOUT, N, M, X, BL, BU, NMAX, EQUATN,     &
                         LINEAR, V, CL, CU, MMAX, EFIRST, LFIRST, NVFRST )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: INPUT, IOUT, N, M, NMAX, MMAX
      LOGICAL :: EFIRST, LFIRST, NVFRST
      REAL ( KIND = wp ) :: X ( NMAX ), BL ( NMAX ), BU ( NMAX )
      REAL ( KIND = wp ) :: V ( MMAX ), CL ( MMAX ), CU ( MMAX )
      LOGICAL :: EQUATN( MMAX ), LINEAR( MMAX )

!  Set up the input data for the constrained optimization tools.

!  Nick Gould, for CGT productions,
!  30th October, 1991.
!  Ingrid Bongartz added option to reorder variables, so that
!  nonlinear variables come first.  Within the nonlinear variables,
!  the smaller set of either the nonlinear objective or nonlinear
!  Jacobian variables appears first.
!  5th November, 1992.

!  Local variables.

      INTEGER :: IALGOR, IPRINT, INFORM, I, IG, J, JG, MEND
      INTEGER :: MEQ, MLIN, NSLACK, NELTYP, NGRTYP
      INTEGER :: II, K, IEL, JWRK, KNDV, NNLIN, NEND
      INTEGER :: ITEMP
      LOGICAL :: FDGRAD, DEBUG, LTEMP
      REAL :: DUM,    CPUTIM
      CHARACTER ( LEN = 8 ) :: PNAME
      CHARACTER ( LEN = 10 ) :: CTEMP
      CHARACTER ( LEN = 10 ) :: CHTEMP
      REAL ( KIND = wp ) :: ATEMP, ZERO
      PARAMETER ( ZERO = 0.0_wp )
      REAL ( KIND = wp ) :: OBFBND( 2 )
      EXTERNAL :: RANGE, CPUTIM
      data%sutime = CPUTIM( DUM )
      data%iout2 = IOUT
      DEBUG = .FALSE.
      DEBUG = DEBUG .AND. IOUT > 0
      IPRINT = 0
      IF ( DEBUG ) IPRINT = 3

!  Input the problem dimensions.

      READ( INPUT, 1001 ) N, data%ng, data%nelnum, data%ngel, data%nvars,      &
         data%nnza, data%ngpvlu, data%nepvlu, NELTYP, NGRTYP
      IF ( N <= 0 ) THEN
         CLOSE( INPUT )
         IF ( IOUT > 0 ) WRITE( IOUT, 2020 )
         STOP
      END IF
      IF ( data%ng <= 0 ) THEN
         CLOSE( INPUT )
         IF ( IOUT > 0 ) WRITE( IOUT, 2030 )
         STOP
      END IF
      IF ( N > NMAX ) THEN
         CLOSE( INPUT )
         IF ( IOUT > 0 ) THEN
            WRITE( IOUT, 2000 ) 'X   ', 'NMAX  ', N - NMAX
            WRITE( IOUT, 2000 ) 'BL  ', 'NMAX  ', N - NMAX
            WRITE( IOUT, 2000 ) 'BU  ', 'NMAX  ', N - NMAX
         END IF
         STOP
      END IF

!  Input the problem type.

      READ( INPUT, 1000 ) IALGOR, PNAME

!  Set useful integer values.

      data%ng1 = data%ng + 1
      data%ngng = data%ng + data%ng
      data%nel1 = data%nelnum + 1

!  Partition the integer workspace.

      ISTADG = 0
      ISTGP = ISTADG + data%ng1
      ISTADA = ISTGP + data%ng1
      ISTAEV = ISTADA + data%ng1
      ISTEP = ISTAEV + data%nel1
      ITYPEG = ISTEP + data%nel1
      KNDOFC = ITYPEG + data%ng
      ITYPEE = KNDOFC + data%ng
      IELING = ITYPEE + data%nelnum
      IELVAR = IELING + data%ngel
      ICNA = IELVAR + data%nvars
      ISTADH = ICNA + data%nnza
      INTVAR = ISTADH + data%nel1
      IVAR = INTVAR + data%nel1
      ICALCF = IVAR + N
      ITYPEV = ICALCF + MAX( data%nelnum, data%ng )
      IWRK = ITYPEV + N
      data%liwork = LIWK - IWRK

!  Ensure there is sufficient room.

      IF ( data%liwork < 0 ) THEN
         CLOSE( INPUT )
         IF ( IOUT > 0 ) WRITE( IOUT, 2000 ) &
             'IWK   ', 'LIWK  ', - data%liwork
         STOP
      END IF

!  Partition the real workspace.

      A = 0
      B = A + data%nnza
      U = B + data%ng
      GPVALU = U + data%ng
      EPVALU = GPVALU + data%ngpvlu
      ESCALE = EPVALU + data%nepvlu
      GSCALE = ESCALE + data%ngel
      VSCALE = GSCALE + data%ng
      GVALS = VSCALE + N
      XT = GVALS + 3 * data%ng
      DGRAD = XT + N
      Q = DGRAD + N
      FT = Q + N
      WRK = FT + data%ng
      data%lwork = LWK - WRK

!  Ensure there is sufficient room.

      IF ( data%lwork < 0 ) THEN
         CLOSE( INPUT )
         IF ( IOUT > 0 ) WRITE( IOUT, 2000 ) &
             'WK   ', 'LWK   ', - data%lwork
         STOP
      END IF

!  Partition the logical workspace.

      INTREP = 0
      GXEQX = INTREP + data%nelnum
      data%lo = GXEQX + data%ngng

!  Ensure there is sufficient room.

      IF ( LLOGIC < data%lo ) THEN
         CLOSE( INPUT )
         IF ( IOUT > 0 ) WRITE( IOUT, 2000 ) &
             'LOGI  ', 'LLOGIC', data%lo - LLOGIC
         STOP
      END IF

!  Partition the character workspace.

      GNAMES = 0
      VNAMES = GNAMES + data%ng
      data%ch = VNAMES + N

!  Ensure there is sufficient room.

      IF ( LCHARA < data%ch + 1 ) THEN
         CLOSE( INPUT )
         IF ( IOUT > 0 ) WRITE( IOUT, 2000 ) &
             'CHA   ', 'LCHARA', data%ch + 1 - LCHARA
         STOP
      END IF

!  Record the lengths of arrays.

      data%lstadg = MAX( 1, data%ng1 )
      data%lstada = MAX( 1, data%ng1 )
      data%lstaev = MAX( 1, data%nel1 )
      data%lkndof = MAX( 1, data%ng )
      data%leling = MAX( 1, data%ngel )
      data%lelvar = MAX( 1, data%nvars )
      data%licna = MAX( 1, data%nnza )
      data%lstadh = MAX( 1, data%nel1 )
      data%lntvar = MAX( 1, data%nel1 )
      data%lcalcf = MAX( 1, data%nelnum, data%ng )
      data%lcalcg = MAX( 1, data%ng )
      data%la = MAX( 1, data%nnza )
      data%lb = MAX( 1, data%ng )
      data%lu = MAX( 1, data%ng )
      data%lescal = MAX( 1, data%ngel )
      data%lgscal = MAX( 1, data%ng )
      data%lvscal = MAX( 1, N )
      data%lft = MAX( 1, data%ng )
      data%lgvals = MAX( 1, data%ng )
      data%lintre = MAX( 1, data%nelnum )
      data%lgxeqx = MAX( 1, data%ngng )
      data%lgpvlu = MAX( 1, data%ngpvlu )
      data%lepvlu = MAX( 1, data%nepvlu )
!     LSTGP = MAX( 1, data%ng1 )
!     LSTEP = MAX( 1, data%nel1 )
!     LTYPEG = MAX( 1, data%ng )
!     LTYPEE = MAX( 1, data%nelnum )
!     LIVAR = MAX( 1, N )
!     LBL = MAX( 1, N )
!     LBU = MAX( 1, N )
!     LX = MAX( 1, N )
!     LXT = MAX( 1, N )
!     LDGRAD = MAX( 1, N )
!     LQ = MAX( 1, N )


!  Print out problem data. input the number of variables, groups,
!  elements and the identity of the objective function group.

      IF ( IALGOR == 2 ) THEN
         READ( INPUT, 1002 ) NSLACK, data%nobjgr
      ELSE
         NSLACK = 0
      END IF
      IF ( DEBUG ) WRITE( IOUT, 1100 ) PNAME, N, data%ng, data%nelnum
      data%pname = PNAME // '  '

!  Input the starting addresses of the elements in each group,
!  of the parameters used for each group and
!  of the nonzeros of the linear element in each group.

      READ( INPUT, 1010 ) ( data%ISTADG( I ), I = 1, data%ng1 )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'ISTADG', &
        ( data%ISTADG( I ), I = 1, data%ng1 )
      READ( INPUT, 1010 ) ( data%ISTGP( I ), I = 1, data%ng1 )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'ISTGP ', &
        ( data%ISTGP( I ), I = 1, data%ng1 )
      READ( INPUT, 1010 ) ( data%ISTADA( I ), I = 1, data%ng1 )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'ISTADA', &
        ( data%ISTADA( I ), I = 1, data%ng1 )

!  Input the starting addresses of the variables and parameters
!  in each element.

      READ( INPUT, 1010 ) ( data%ISTAEV( I ), I = 1, data%nel1 )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'ISTAEV', &
        ( data%ISTAEV( I ), I = 1, data%nel1 )
      READ( INPUT, 1010 ) ( data%ISTEP( I ), I = 1, data%nel1 )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'ISTEP ', &
        ( data%ISTEP( I ), I = 1, data%nel1 )

!  Input the group type of each group

      READ( INPUT, 1010 ) ( data%ITYPEG( I ), I = 1, data%ng )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'ITYPEG', &
        ( data%ITYPEG( I ), I = 1, data%ng )
      IF ( IALGOR >= 2 ) THEN
         READ( INPUT, 1010 ) ( data%KNDOFC( I ), I = 1, data%ng )
         IF ( DEBUG ) WRITE( IOUT, 1110 ) 'KNDOFC', &
        ( data%KNDOFC( I ), I = 1, data%ng )
      END IF

!  Input the element type of each element

      READ( INPUT, 1010 ) ( data%ITYPEE( I ), I = 1, data%nelnum )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'ITYPEE', &
        ( data%ITYPEE( I ), I = 1, data%nelnum )

!  Input the number of internal variables for each element.

      READ( INPUT, 1010 ) ( data%INTVAR( I ), I = 1, data%nelnum )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'INTVAR', &
        ( data%INTVAR( I ), I = 1, data%nelnum )

!  Input the identity of each individual element.

      READ( INPUT, 1010 ) ( data%IELING( I ), I = 1, data%ngel )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'IELING', &
        ( data%IELING( I ), I = 1, data%ngel )

!  Input the variables in each group's elements.

      data%nvars = data%ISTAEV( data%nel1 ) - 1
      READ( INPUT, 1010 ) ( data%IELVAR( I ), I = 1, data%nvars )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'IELVAR', &
        ( data%IELVAR( I ), I = 1, data%nvars )

!  Input the column addresses of the nonzeros in each linear element.

      READ( INPUT, 1010 ) ( data%ICNA( I ), I = 1, data%nnza )
      IF ( DEBUG ) WRITE( IOUT, 1110 ) 'ICNA  ', &
        ( data%ICNA( I ), I = 1, data%nnza )

!  Input the values of the nonzeros in each linear element, the
!  constant term in each group, the lower and upper bounds on
!  the variables.

      READ( INPUT, 1020 ) ( data%A( I ), I = 1, data%nnza )
      IF ( DEBUG ) WRITE( IOUT, 1120 ) 'A     ', &
        ( data%A( I ), I = 1, data%nnza )
      READ( INPUT, 1020 ) ( data%B( I ), I = 1, data%ng )
      IF ( DEBUG ) WRITE( IOUT, 1120 ) 'B     ', &
        ( data%B( I ), I = 1, data%ng )
      IF ( IALGOR <= 2 ) THEN
         READ( INPUT, 1020 ) ( BL( I ), I = 1, N )
         IF ( DEBUG ) WRITE( IOUT, 1120 ) 'BL    ', ( BL( I ), I = 1, N )
         READ( INPUT, 1020 ) ( BU( I ), I = 1, N )
         IF ( DEBUG ) WRITE( IOUT, 1120 ) 'BU    ', ( BU( I ), I = 1, N )
      ELSE

!  Use GVALS and FT as temporary storage for the constraint bounds.

         READ( INPUT, 1020 ) ( BL( I ), I = 1, N ), &
           ( data%GVALS( I ), I = 1, data%ng )
         IF ( DEBUG ) WRITE( IOUT, 1120 ) 'BL    ', &
           ( BL( I ), I = 1, N ), ( data%GVALS( I ), I = 1, data%ng )
         READ( INPUT, 1020 ) ( BU( I ), I = 1, N ), &
           ( data%FT( I ), I = 1, data%ng )
         IF ( DEBUG ) WRITE( IOUT, 1120 ) 'BU    ', &
           ( BU( I ), I = 1, N ), ( data%FT( I ), I = 1, data%ng )
      END IF

!   Input the starting point for the minimization.

      READ( INPUT, 1020 ) ( X( I ), I = 1, N )
      IF ( DEBUG ) WRITE( IOUT, 1120 ) 'X     ', ( X( I ), I = 1, N )
      IF ( IALGOR >= 2 ) THEN
         READ( INPUT, 1020 )( data%U( I ), I = 1, data%ng )
         IF ( DEBUG ) WRITE( IOUT, 1120 ) 'U     ', &
            ( data%U( I ), I = 1, data%ng )
      END IF

!  Input the parameters in each group.

      READ( INPUT, 1020 ) ( data%GPVALU( I ), I = 1, data%ngpvlu )
      IF ( DEBUG ) WRITE( IOUT, 1120 ) 'GPVALU', &
        ( data%GPVALU( I ), I = 1, data%ngpvlu )

!  Input the parameters in each individual element.

      READ( INPUT, 1020 ) ( data%EPVALU( I ), I = 1, data%nepvlu )
      IF ( DEBUG ) WRITE( IOUT, 1120 ) 'EPVALU', &
        ( data%EPVALU( I ), I = 1, data%nepvlu )

!  Input the scale factors for the nonlinear elements.

      READ( INPUT, 1020 ) ( data%ESCALE( I ), I = 1, data%ngel )
      IF ( DEBUG ) WRITE( IOUT, 1120 ) 'ESCALE', &
        ( data%ESCALE( I ), I = 1, data%ngel )

!  Input the scale factors for the groups.

      READ( INPUT, 1020 ) ( data%GSCALE( I ), I = 1, data%ng )
      IF ( DEBUG ) WRITE( IOUT, 1120 ) 'GSCALE', &
        ( data%GSCALE( I ), I = 1, data%ng )

!  Input the scale factors for the variables.

      READ( INPUT, 1020 ) ( data%VSCALE( I ), I = 1, N )
      IF ( DEBUG ) WRITE( IOUT, 1120 ) 'VSCALE', &
        ( data%VSCALE( I ), I = 1, N )

!  Input the lower and upper bounds on the objective function.

      READ( INPUT, 1080 ) OBFBND( 1 ), OBFBND( 2 )
      IF ( DEBUG ) WRITE( IOUT, 1180 ) 'OBFBND', &
          OBFBND( 1 ), OBFBND( 2 )

!  Input a logical array which says whether an element has internal
!  variables.

      READ( INPUT, 1030 ) ( data%INTREP( I ), I = 1, data%nelnum )
      IF ( DEBUG ) WRITE( IOUT, 1130 ) 'INTREP', &
        ( data%INTREP( I ), I = 1, data%nelnum )

!  Input a logical array which says whether a group is trivial.

      READ( INPUT, 1030 ) ( data%GXEQX( I ), I = 1, data%ng )
      IF ( DEBUG ) WRITE( IOUT, 1130 ) 'GXEQX ', &
        ( data%GXEQX( I ), I = 1, data%ng )

!  Input the names given to the groups and to the variables.

      READ( INPUT, 1040 ) ( data%GNAMES( I ), I = 1, data%ng )
      IF ( DEBUG ) WRITE( IOUT, 1140 ) 'GNAMES', &
        ( data%GNAMES( I ), I = 1, data%ng )
      READ( INPUT, 1040 ) ( data%VNAMES( I ), I = 1, N )
      IF ( DEBUG ) WRITE( IOUT, 1140 ) 'VNAMES', &
        ( data%VNAMES( I ), I = 1, N )

!  Dummy input for the names given to the element and group types.

      READ( INPUT, 1040 ) ( CHTEMP, I = 1, NELTYP )
      READ( INPUT, 1040 ) ( CHTEMP, I = 1, NGRTYP )

!  Input the type of each variable.

      READ( INPUT, 1010 ) ( data%ITYPEV( I ), I = 1, N )

!  Consider which groups are constraints. Of these, decide which are
!  equations, which are linear, allocate starting values for the
!  Lagrange multipliers and set lower and upper bounds on any
!  inequality constraints. Reset KNDOFC to point to the list of
!  constraint groups.

      M = 0
      DO 10 I = 1, data%ng
         IF ( data%KNDOFC( I ) == 1 ) THEN
            data%KNDOFC( I ) = 0
         ELSE
            M = M + 1
            IF ( M <= MMAX ) THEN
               V ( M ) = data%U( I )
               LINEAR( M ) = data%GXEQX( I ) .AND. data%ISTADG( I ) &
                             >= data%ISTADG( I + 1 )
               IF ( data%KNDOFC( I ) == 2 ) THEN
                  EQUATN( M ) = .TRUE.
                  CL ( M ) = ZERO
                  CU ( M ) = ZERO
               ELSE
                  EQUATN( M ) = .FALSE.
                  CL ( M ) = data%GVALS( I )
                  CU ( M ) = data%FT( I )
               END IF
            END IF
            data%KNDOFC( I ) = M
         END IF
   10 CONTINUE
      IF ( M == 0 .AND. IOUT > 0 ) WRITE( IOUT, 2010 )
      IF ( M > MMAX ) THEN
         CLOSE( INPUT )
         IF ( IOUT > 0 ) THEN
            WRITE( IOUT, 2000 ) 'V     ', 'MMAX  ', M - MMAX
            WRITE( IOUT, 2000 ) 'CL    ', 'MMAX  ', M - MMAX
            WRITE( IOUT, 2000 ) 'CU    ', 'MMAX  ', M - MMAX
            WRITE( IOUT, 2000 ) 'EQUATN', 'MMAX  ', M - MMAX
            WRITE( IOUT, 2000 ) 'LINEAR', 'MMAX  ', M - MMAX
         END IF
         STOP
      END IF
!                            the number of variables and constraints, resp.
      data%numvar = N
      data%numcon = M
      IF ( NVFRST ) THEN

!  Ensure there is sufficient room in IWK to reorder variables.

         IF ( data%liwork < 2*N ) THEN
            CLOSE( INPUT )
            IF ( IOUT > 0 ) WRITE( IOUT, 2000 ) &
                'IWK   ', 'LIWK  ', data%liwork - 2*N
            STOP
         END IF
         KNDV = IWRK + 1 
         JWRK = KNDV + N

!  Initialize JWRK and KNDV.

         DO 20 J = 1, N
            data%IWORK( KNDV + J ) = 0
            data%IWORK( JWRK + J ) = J
   20    CONTINUE

!  Now identify and count nonlinear variables.
!  Keep separate counts for nonlinear objective and Jacobian variables.
!  data%IWORK( KNDV + J ) = 0 ==> J linear everywhere
!  data%IWORK( KNDV + J ) = 1 ==> J linear in objective, nonlinear in constraints
!  data%IWORK( KNDV + J ) = 2 ==> J linear in constraints, nonlinear in objective
!  data%IWORK( KNDV + J ) = 3 ==> J nonlinear everywhere

         NNLIN = 0
         data%nnov = 0
         data%nnjv = 0
         DO 60 IG = 1, data%ng
            I = data%KNDOFC( IG )
            DO 40 II = data%ISTADG( IG ), data%ISTADG( IG + 1 ) - 1
               IEL = data%IELING( II )
               DO 30 K = data%ISTAEV( IEL ), data%ISTAEV( IEL + 1) - 1
                  J = data%IELVAR( K )
                  IF ( I > 0 ) THEN
                     IF ( data%IWORK( KNDV + J ) == 0 ) THEN
                        data%IWORK( KNDV + J ) = 1
                        data%nnjv = data%nnjv + 1
                        NNLIN = NNLIN + 1
                     ELSE IF ( data%IWORK( KNDV + J ) == 2 ) THEN
                        data%IWORK( KNDV + J ) = 3
                        data%nnjv = data%nnjv + 1
                     END IF
                  ELSE
                     IF ( data%IWORK( KNDV + J ) == 0 ) THEN
                        data%IWORK( KNDV + J ) = 2
                        data%nnov = data%nnov + 1
                        NNLIN = NNLIN + 1
                     ELSE IF ( data%IWORK( KNDV + J ) == 1 ) THEN
                        data%IWORK( KNDV + J ) = 3
                        data%nnov = data%nnov + 1
                     END IF
                  END IF
   30          CONTINUE
   40       CONTINUE
            IF ( .NOT. data%GXEQX( IG ) ) THEN
               DO 50 II = data%ISTADA( IG ), data%ISTADA( IG + 1 ) - 1
                  J = data%ICNA( II )
                  IF ( I > 0 ) THEN
                     IF ( data%IWORK( KNDV + J ) == 0 ) THEN
                        data%IWORK( KNDV + J ) = 1
                        data%nnjv = data%nnjv + 1
                        NNLIN = NNLIN + 1
                     ELSE IF ( data%IWORK( KNDV + J ) == 2 ) THEN
                        data%IWORK( KNDV + J ) = 3
                        data%nnjv = data%nnjv + 1
                     END IF
                  ELSE
                     IF ( data%IWORK( KNDV + J ) == 0 ) THEN
                        data%IWORK( KNDV + J ) = 2
                        data%nnov = data%nnov + 1
                        NNLIN = NNLIN + 1
                     ELSE IF ( data%IWORK( KNDV + J ) == 1 ) THEN
                        data%IWORK( KNDV + J ) = 3
                        data%nnov = data%nnov + 1
                     END IF
                  END IF
   50          CONTINUE
            END IF
   60    CONTINUE
         IF ( NNLIN == 0 .OR. ( data%nnov == N .AND. data%nnjv == N ) ) &
            GO TO 600
         IF ( NNLIN == N ) GO TO 500

!  Reorder the variables so that all nonlinear variables occur before
!  the linear ones.

         NEND = N

!  Run forward through the variables until a linear variable
!  is encountered.

         DO 420 I = 1, N
            IF ( I > NEND ) GO TO 430
            IF ( data%IWORK( KNDV + I ) == 0 ) THEN

!  Variable I is linear. Now, run backwards through the
!  variables until a nonlinear one is encountered.

               DO 410 J = NEND, I, - 1
                  IF ( data%IWORK( KNDV + J ) > 0 ) THEN 
                     NEND = J - 1

!  Interchange the data for variables I and J.

                     ITEMP = data%IWORK( JWRK + I )
                     data%IWORK( JWRK + I ) = data%IWORK( JWRK + J )
                     data%IWORK( JWRK + J ) = ITEMP
                     ITEMP = data%IWORK( KNDV + I )
                     data%IWORK( KNDV + I ) = data%IWORK( KNDV + J )
                     data%IWORK( KNDV + J ) = ITEMP
                     ATEMP = BL ( I )
                     BL ( I ) = BL ( J )
                     BL ( J ) = ATEMP
                     ATEMP = BU ( I )
                     BU ( I ) = BU ( J )
                     BU ( J ) = ATEMP
                     ATEMP = X ( I )
                     X ( I ) = X ( J )
                     X ( J ) = ATEMP
                     ATEMP = data%VSCALE( I )
                     data%VSCALE( I ) = data%VSCALE( J )
                     data%VSCALE( J ) = ATEMP
                     CTEMP = data%VNAMES( I )
                     data%VNAMES( I ) = data%VNAMES( J )
                     data%VNAMES( J ) = CTEMP
                     GO TO 420
                  END IF
  410          CONTINUE
               GO TO 430
            END IF
  420    CONTINUE
  430    CONTINUE 

!  Change entries in IELVAR and ICNA to reflect reordering of variables.

         DO 440 I = 1, data%nvars
            J = data%IELVAR( I )
            data%IELVAR( I ) = data%IWORK( JWRK + J ) 
  440    CONTINUE
         DO 450 I = 1, data%nnza
            J = data%ICNA( I )
            data%ICNA( I ) = data%IWORK( JWRK + J )
  450    CONTINUE
         DO 460 J = 1, N
            data%IWORK( JWRK + J ) = J
  460    CONTINUE
  500    CONTINUE
         IF ( ( data%nnov == NNLIN .AND. data%nnjv == NNLIN )  &
            .OR. ( data%nnov == 0 ) .OR. ( data%nnjv == 0 ) ) GO TO 600

!  Reorder the nonlinear variables so that the smaller set (nonlinear
!  objective or nonlinear Jacobian) occurs at the beginning of the 
!  larger set.

         NEND = NNLIN
         IF ( data%nnjv <= data%nnov ) THEN

!  Put the nonlinear Jacobian variables first.
!  Reset data%nnov to indicate all nonlinear variables are treated as
!  nonlinear objective variables.

            data%nnov = NNLIN
            DO 520 I = 1, NNLIN 
               IF ( I > NEND ) GO TO 530
               IF ( data%IWORK( KNDV + I ) == 2 ) THEN

!  Variable I is linear in the Jacobian. Now, run backwards through the 
!  variables until a nonlinear Jacobian variable is encountered.

                  DO 510 J = NEND, I, - 1
                     IF ( data%IWORK( KNDV + J ) == 1  &
                        .OR. data%IWORK( KNDV + J ) == 3 ) THEN 
                        NEND = J - 1

!  Interchange the data for variables I and J.

                        ITEMP = data%IWORK( JWRK + I )
                        data%IWORK( JWRK + I ) = data%IWORK( JWRK + J )
                        data%IWORK( JWRK + J ) = ITEMP
                        ITEMP = data%IWORK( KNDV + I )
                        data%IWORK( KNDV + I ) = data%IWORK( KNDV + J )
                        data%IWORK( KNDV + J ) = ITEMP
                        ATEMP = BL ( I )
                        BL ( I ) = BL ( J )
                        BL ( J ) = ATEMP
                        ATEMP = BU ( I )
                        BU ( I ) = BU ( J )
                        BU ( J ) = ATEMP
                        ATEMP = X ( I )
                        X ( I ) = X ( J )
                        X ( J ) = ATEMP
                        ATEMP = data%VSCALE( I )
                        data%VSCALE( I ) = data%VSCALE( J )
                        data%VSCALE( J ) = ATEMP
                        CTEMP = data%VNAMES( I )
                        data%VNAMES( I ) = data%VNAMES( J )
                        data%VNAMES( J ) = CTEMP
                        GO TO 520
                     END IF
  510             CONTINUE
                  GO TO 530
               END IF
  520       CONTINUE
  530       CONTINUE
         ELSE

!  Put the nonlinear objective variables first.
!  Reset data%nnjv to indicate all nonlinear variables are treated as
!  nonlinear Jacobian variables.

            data%nnjv = NNLIN
            DO 550 I = 1, NNLIN 
               IF ( I > NEND ) GO TO 560
               IF ( data%IWORK( KNDV + I ) == 1 ) THEN

!  Variable I is linear in the objective. Now, run backwards through the 
!  variables until a nonlinear objective variable is encountered.

                  DO 540 J = NEND, I, - 1
                     IF ( data%IWORK( KNDV + J ) > 1 ) THEN 
                        NEND = J - 1

!  Interchange the data for variables I and J.

                        ITEMP = data%IWORK( JWRK + I )
                        data%IWORK( JWRK + I ) = data%IWORK( JWRK + J )
                        data%IWORK( JWRK + J ) = ITEMP
                        ITEMP = data%IWORK( KNDV + I )
                        data%IWORK( KNDV + I ) = data%IWORK( KNDV + J )
                        data%IWORK( KNDV + J ) = ITEMP
                        ATEMP = BL ( I )
                        BL ( I ) = BL ( J )
                        BL ( J ) = ATEMP
                        ATEMP = BU ( I )
                        BU ( I ) = BU ( J )
                        BU ( J ) = ATEMP
                        ATEMP = X ( I )
                        X ( I ) = X ( J )
                        X ( J ) = ATEMP
                        ATEMP = data%VSCALE( I )
                        data%VSCALE( I ) = data%VSCALE( J )
                        data%VSCALE( J ) = ATEMP
                        CTEMP = data%VNAMES( I )
                        data%VNAMES( I ) = data%VNAMES( J )
                        data%VNAMES( J ) = CTEMP
                        GO TO 550
                     END IF
  540             CONTINUE
                  GO TO 560
               END IF
  550       CONTINUE
  560       CONTINUE
         END IF

!  Change entries in IELVAR and ICNA to reflect reordering of variables.

         DO 580 I = 1, data%nvars
            J = data%IELVAR( I )
            data%IELVAR( I ) = data%IWORK( JWRK + J ) 
  580    CONTINUE
         DO 590 I = 1, data%nnza
            J = data%ICNA( I )
            data%ICNA( I ) = data%IWORK( JWRK + J )
  590    CONTINUE
  600    CONTINUE
      END IF

!  Partition the workspace arrays data%FUVALS, IWK and WK. Initialize
!  certain portions of IWK.

      data%firstg = .TRUE.
      FDGRAD = .FALSE.
      CALL DINITW( N, data%ng, data%nelnum, data%IELING, data%leling,          &
          data%ISTADG, data%lstadg, data%IELVAR, data%lelvar, data%ISTAEV,     &
          data%lstaev, data%INTVAR, data%lntvar, data%ISTADH, data%lstadh,     &
          data%ICNA, data%licna, data%ISTADA, data%lstada, data%ITYPEE,        &
          data%lintre, data%GXEQX, data%lgxeqx, data%INTREP, data%lintre,      &
          LFUVAL, data%altriv, .TRUE., FDGRAD, data%lfxi,LGXI,LHXI,LGGFX,      &
          data%ldx, data%lgrjac, data%lqgrad, data%lbreak, data%lp, data%lxcp, &
          data%lx0, data%lgx0, data%ldeltx, data%lbnd, data%lwkstr,            &
          data%lsptrs, data%lselts, data%lindex, data%lswksp, data%lstagv,     &
          data%lstajc, data%liused, data%lfreec, data%lnnonz, data%lnonz2,     &
          data%lsymmd, data%lsymmh, data%lslgrp, data%lsvgrp, data%lgcolj,     &
          data%lvaljr, data%lsend,  data%lnptrs, data%lnelts, data%lnndex,     &
          data%lnwksp, data%lnstgv, data%lnstjc, data%lniuse,  data%lnfrec,    &
          data%lnnnon, data%lnnno2, data%lnsymd, data%lnsymh, data%lnlgrp,     &
          data%lnvgrp, data%lngclj, data%lnvljr, data%lnqgrd, data%lnbrak,     &
          data%lnp, data%lnbnd, data%lnfxi,  data%lngxi,  data%lnguvl,         &
          data%lnhxi,  data%lnhuvl, data%lnggfx, data%lndx, data%lngrjc,       &
          data%liwk2, data%lwk2, data%maxsin, data%ninvar, data%ntype,         &
          data%nsets, data%maxsel, data%lstype, data%lsswtr, data%lssiwt,      &
          data%lsiwtr, data%lswtra, data%lntype, data%lnswtr, data%lnsiwt,     &
          data%lniwtr, data%lnwtra, data%lsiset, data%lssvse, data%lniset,     &
          data%lnsvse, RANGE, data%IWORK(IWRK + 1), LIWORK, data%WRK, LWORK,   &
          IPRINT, IOUT, INFORM )
      IF ( INFORM /= 0 ) STOP

!  Shift the starting addresses for the real workspace relative to WRK.

      data%lqgrad = data%lqgrad + WRK
      data%lbreak = data%lbreak + WRK
      data%lp = data%lp + WRK
      data%lxcp = data%lxcp + WRK
      data%lx0 = data%lx0 + WRK
      data%lgx0 = data%lgx0 + WRK
      data%ldeltx = data%ldeltx + WRK
      data%lbnd = data%lbnd + WRK
      data%lswtra = data%lswtra + WRK
      data%lwkstr = data%lwkstr + WRK

!  Shift the starting addresses for the integer workspace relative
!  to IWRK.

      data%lsptrs = data%lsptrs + IWRK
      data%lselts = data%lselts + IWRK
      data%lindex = data%lindex + IWRK
      data%lswksp = data%lswksp + IWRK
      data%lstagv = data%lstagv + IWRK
      data%lstajc = data%lstajc + IWRK
      data%liused = data%liused + IWRK
      data%lfreec = data%lfreec + IWRK
      data%lnnonz = data%lnnonz + IWRK
      data%lnonz2 = data%lnonz2 + IWRK
      data%lsymmd = data%lsymmd + IWRK
      data%lsymmh = data%lsymmh + IWRK
      data%lslgrp = data%lslgrp + IWRK
      data%lsvgrp = data%lsvgrp + IWRK
      data%lgcolj = data%lgcolj + IWRK
      data%lvaljr = data%lvaljr + IWRK
      data%lstype = data%lstype + IWRK
      data%lsswtr = data%lsswtr + IWRK
      data%lssiwt = data%lssiwt + IWRK
      data%lsiwtr = data%lsiwtr + IWRK
      data%lsiset = data%lsiset + IWRK
      data%lssvse = data%lssvse + IWRK
      data%lsend = data%lsend + IWRK
      IF ( .NOT. ( EFIRST .OR. LFIRST ) ) GOTO 340
!                            to RETURN if there are no constraints.
      IF ( M == 0 ) GOTO 340

!  Record which group is associated with each constraint.

      IF ( M > data%liwk2 ) THEN
         WRITE( IOUT, 2040 )
         STOP
      END IF
      MEQ = 0
      MLIN = 0
      DO 100 IG = 1, data%ng
         I = data%KNDOFC( IG )
         IF ( I > 0 ) THEN
            data%IWORK( data%lsend + I ) = IG
            IF ( EQUATN( I ) ) MEQ = MEQ + 1
            IF ( LINEAR( I ) ) MLIN = MLIN + 1
         END IF
  100 CONTINUE
      IF ( LFIRST ) THEN
         IF ( MLIN == 0 .OR. MLIN == M ) GO TO 130

!  Reorder the constraints so that the linear constraints occur before the
!  nonlinear ones.

         MEND = M

!  Run forward through the constraints until a nonlinear constraint
!  is encountered.

         DO 120 I = 1, M
            IF ( I > MEND ) GO TO 130
            IG = data%IWORK( data%lsend + I )
!              write(6,*) ' group ', IG, ' type ', I, ' equal? ',
!     *                     EQUATN( I )
            IF ( .NOT. LINEAR( I ) ) THEN

!  Constraint I is nonlinear. Now, run backwards through the
!  constraints until a linear one is encountered.

               DO 110 J = MEND, I, - 1
                  JG = data%IWORK( data%lsend + J )
!                 write(6,*) ' group ', JG, ' type ', J,
!     *                      ' linear? ', LINEAR( J )
                  IF ( LINEAR( J ) ) THEN
!                    write(6,*) ' swaping constraints ', I,
!     *                         ' and ', J
                     MEND = J - 1

!  Interchange the data for constraints I and J.

                     data%IWORK( data%lsend + I ) = JG
                     data%IWORK( data%lsend + J ) = IG
                     data%KNDOFC( IG ) = J
                     data%KNDOFC( JG ) = I
                     LTEMP = LINEAR( I )
                     LINEAR( I ) = LINEAR( J )
                     LINEAR( J ) = LTEMP
                     LTEMP = EQUATN( I )
                     EQUATN( I ) = EQUATN( J )
                     EQUATN( J ) = LTEMP
                     ATEMP = V ( I )
                     V ( I ) = V ( J )
                     V ( J ) = ATEMP
                     ATEMP = CL ( I )
                     CL ( I ) = CL ( J )
                     CL ( J ) = ATEMP
                     ATEMP = CU ( I )
                     CU ( I ) = CU ( J )
                     CU ( J ) = ATEMP
                     GO TO 120
                  END IF
  110          CONTINUE
               GO TO 130
            END IF
  120    CONTINUE
  130    CONTINUE
         IF ( EFIRST ) THEN
            IF ( MEQ == 0 .OR. MEQ == M ) GO TO 260

!  Reorder the linear constraints so that the equations occur before
!  the inequalities.

            MEND = MLIN
            DO 220 I = 1, MLIN
               IF ( I > MEND ) GO TO 230
               IG = data%IWORK( data%lsend + I )
!                 write(6,*) ' group ', IG, ' type ', I, ' equation? ',
!     *                        EQUATN( I )
               IF ( .NOT. EQUATN( I ) ) THEN

!  Constraint I is an inequality. Now, run backwards through the
!  constraints until an equation is encountered.

                  DO 210 J = MEND, I, - 1
                     JG = data%IWORK( data%lsend + J )
!                    write(6,*) ' group ', JG, ' type ', J,
!     *                         ' equation? ', EQUATN( J )
                     IF ( EQUATN( J ) ) THEN
!                       write(6,*) ' swaping constraints ', I,
!     *                            ' and ', J
                        MEND = J - 1

!  Interchange the data for constraints I and J.

                        data%IWORK( data%lsend + I ) = JG
                        data%IWORK( data%lsend + J ) = IG
                        data%KNDOFC( IG ) = J
                        data%KNDOFC( JG ) = I
                        LTEMP = LINEAR( I )
                        LINEAR( I ) = LINEAR( J )
                        LINEAR( J ) = LTEMP
                        LTEMP = EQUATN( I )
                        EQUATN( I ) = EQUATN( J )
                        EQUATN( J ) = LTEMP
                        ATEMP = V ( I )
                        V ( I ) = V ( J )
                        V ( J ) = ATEMP
                        ATEMP = CL ( I )
                        CL ( I ) = CL ( J )
                        CL ( J ) = ATEMP
                        ATEMP = CU ( I )
                        CU ( I ) = CU ( J )
                        CU ( J ) = ATEMP
                        GO TO 220
                     END IF
  210             CONTINUE
                  GO TO 230
               END IF
  220       CONTINUE
  230       CONTINUE

!  Reorder the nonlinear constraints so that the equations occur
!  before the inequalities.

            MEND = M
            DO 250 I = MLIN + 1, M
               IF ( I > MEND ) GO TO 260
               IG = data%IWORK( data%lsend + I )
!                 write(6,*) ' group ', IG, ' type ', I, ' equation? ',
!     *                        EQUATN( I )
               IF ( .NOT. EQUATN( I ) ) THEN

!  Constraint I is an inequality. Now, run backwards through the
!  constraints until an equation is encountered.

                  DO 240 J = MEND, I, - 1
                     JG = data%IWORK( data%lsend + J )
!                    write(6,*) ' group ', JG, ' type ', J, &
!                             ' equation? ', EQUATN( J )
                     IF ( EQUATN( J ) ) THEN
!                       write(6,*) ' swaping constraints ', I, ' and ', J
                        MEND = J - 1

!  Interchange the data for constraints I and J.

                        data%IWORK( data%lsend + I ) = JG
                        data%IWORK( data%lsend + J ) = IG
                        data%KNDOFC( IG ) = J
                        data%KNDOFC( JG ) = I
                        LTEMP = LINEAR( I )
                        LINEAR( I ) = LINEAR( J )
                        LINEAR( J ) = LTEMP
                        LTEMP = EQUATN( I )
                        EQUATN( I ) = EQUATN( J )
                        EQUATN( J ) = LTEMP
                        ATEMP = V ( I )
                        V ( I ) = V ( J )
                        V ( J ) = ATEMP
                        ATEMP = CL ( I )
                        CL ( I ) = CL ( J )
                        CL ( J ) = ATEMP
                        ATEMP = CU ( I )
                        CU ( I ) = CU ( J )
                        CU ( J ) = ATEMP
                        GO TO 250
                     END IF
  240             CONTINUE
                  GO TO 260
               END IF
  250       CONTINUE
  260       CONTINUE
         END IF
      ELSE
         IF ( EFIRST ) THEN
            IF ( MEQ == 0 .OR. MEQ == M ) GO TO 330

!  Reorder the constraints so that the equations occur before the
!  inequalities.

            MEND = M
            DO 320 I = 1, M
               IF ( I > MEND ) GO TO 330
               IG = data%IWORK( data%lsend + I )
!                 write(6,*) ' group ', IG, ' type ', I, ' equation? ',
!     *                        EQUATN( I )
               IF ( .NOT. EQUATN( I ) ) THEN

!  Constraint I is an inequality. Now, run backwards through the
!  constraints until an equation is encountered.

                  DO 310 J = MEND, I, - 1
                     JG = data%IWORK( data%lsend + J )
!                    write(6,*) ' group ', JG, ' type ', J,                    &
!                               ' equation? ', EQUATN( J )
                     IF ( EQUATN( J ) ) THEN
!                       write(6,*) ' swaping constraints ', I,' and ', J
                        MEND = J - 1

!  Interchange the data for constraints I and J.

                        data%IWORK( data%lsend + I ) = JG
                        data%IWORK( data%lsend + J ) = IG
                        data%KNDOFC( IG ) = J
                        data%KNDOFC( JG ) = I
                        LTEMP = LINEAR( I )
                        LINEAR( I ) = LINEAR( J )
                        LINEAR( J ) = LTEMP
                        LTEMP = EQUATN( I )
                        EQUATN( I ) = EQUATN( J )
                        EQUATN( J ) = LTEMP
                        ATEMP = V ( I )
                        V ( I ) = V ( J )
                        V ( J ) = ATEMP
                        ATEMP = CL ( I )
                        CL ( I ) = CL ( J )
                        CL ( J ) = ATEMP
                        ATEMP = CU ( I )
                        CU ( I ) = CU ( J )
                        CU ( J ) = ATEMP
                        GO TO 320
                     END IF
  310             CONTINUE
                  GO TO 330
               END IF
  320       CONTINUE
  330       CONTINUE
         END IF
      END IF

!  Initialize the performance counters and variables

 340  CONTINUE
      data%nc2of = 0
      data%nc2og = 0
      data%nc2oh = 0
      data%nc2cf = 0
      data%nc2cg = 0
      data%nc2ch = 0
      data%nhvpr = 0
      data%pnc = M
      data%sttime = CPUTIM( DUM )
      data%sutime = data%sttime - data%sutime

      RETURN

!  Non-executable statements.

 1000 FORMAT( I2, A8 )
 1001 FORMAT( 10I8 )
 1002 FORMAT( 2I8 )
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
 2000 FORMAT( /, ' ** SUBROUTINE CSETUP: array length ', A6, &
              ' too small.', /, ' -- Miminimization abandoned.', &
              /, ' -- Increase the parameter ', A6, ' by at least ', I8, &
                 ' and restart.' )
 2010 FORMAT( /, ' ** SUBROUTINE CSETUP: ** Warning. The problem has', &
                 ' no general constraints. ', /, &
                 ' Other tools may be preferable' )
 2020 FORMAT( /, ' ** SUBROUTINE CSETUP: the problem uses no variables.' &
, ' Execution terminating ' )
 2030 FORMAT( /, ' ** SUBROUTINE CSETUP: the problem is vacuous.', &
                 ' Execution terminating ' )
 2040 FORMAT( ' ** SUBROUTINE CSETUP: Increase the size of IWK ' )

!  End of CSETUP.

      END
