! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE USETUP( data, INPUT, IOUT, N, X, BL, BU, NMAX )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: INPUT, IOUT, N, NMAX
      REAL ( KIND = wp ) :: X( NMAX ), BL( NMAX ), BU( NMAX )

!  Set up the input data for the remaining unconstrained
!  optimization tools.

!  Nick Gould, for CGT productions,
!  30th October, 1991.

!  local variables.

      INTEGER :: IALGOR, IPRINT, INFORM, I, NSLACK
      LOGICAL :: FDGRAD, DEBUG
      REAL :: DUM,    CPUTIM
      REAL ( KIND = wp ) :: OBFBND( 2 )
      CHARACTER ( LEN = 8 ) :: PNAME
      CHARACTER ( LEN = 10 ) :: CHTEMP
      EXTERNAL :: RANGE, CPUTIM
      data%sutime = CPUTIM( DUM )
      data%iout2 = IOUT
      DEBUG = .FALSE.
      DEBUG = DEBUG .AND. IOUT > 0
      IPRINT = 0
      IF ( DEBUG ) IPRINT = 3

!  Input the problem dimensions.

      READ( INPUT, 1001 ) N, data%ng, data%nelnum, data%ngel, data%nvars, data%nnza, data%ngpvlu, &
                          data%nepvlu, NELTYP, NGRTYP
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
         READ( INPUT, 1010 )( data%KNDOFC( I ), I = 1, data%ng )
         IF ( DEBUG ) WRITE( IOUT, 1110 ) 'KNDOFC', &
 ( data%KNDOFC( I ), I = 1, data%ng )
         DO 10 I = 1, data%ng
            IF ( ABS( data%KNDOFC( I ) ) >= 2 ) THEN
               CLOSE( INPUT )
               IF ( IOUT > 0 ) WRITE( IOUT, 2010 )
               STOP
            END IF
   10    CONTINUE
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
!  the variables and the starting point for the minimization.

      READ( INPUT, 1020 ) ( data%A( I ), I = 1, data%nnza )
      IF ( DEBUG ) WRITE( IOUT, 1120 ) 'A     ', &
 ( data%A( I ), I = 1, data%nnza )
      READ( INPUT, 1020 ) ( data%B( I ), I = 1, data%ng )
      IF ( DEBUG ) WRITE( IOUT, 1120 ) 'B     ', &
 ( data%B( I ), I = 1, data%ng )
      IF ( IALGOR <= 2 ) THEN
         READ( INPUT, 1020 ) ( BL( I ), I = 1, N )
         IF ( DEBUG ) WRITE( IOUT, 1120 ) 'BL    ', &
 ( BL( I ), I = 1, N )
         READ( INPUT, 1020 ) ( BU( I ), I = 1, N )
         IF ( DEBUG ) WRITE( IOUT, 1120 ) 'BU    ', &
 ( BU( I ), I = 1, N )
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
      READ( INPUT, 1020 ) ( X( I ), I = 1, N )
      IF ( DEBUG ) WRITE( IOUT, 1120 ) 'X     ', &
 ( X( I ), I = 1, N )
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
!  varaiables.

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
      CLOSE( INPUT )

      data%numvar = N


!  Partition the workspace arrays data%FUVALS, IWK and WK. Initialize
!  certain portions of IWK.

      data%firstg = .TRUE.
      FDGRAD = .FALSE.
      CALL DINITW( N, data%ng, data%nelnum, data%IELING, data%leling, data%ISTADG, &
          data%lstadg, data%IELVAR, data%lelvar, data%ISTAEV, data%lstaev, &
          data%INTVAR, data%lntvar, data%ISTADH, data%lstadh, &
          data%ICNA, data%licna, data%ISTADA, data%lstada, &
          data%ITYPEE, data%lintre, &
          data%GXEQX, data%lgxeqx, data%INTREP, data%lintre, &
          LFUVAL, data%altriv, .TRUE., FDGRAD, data%lfxi,LGXI,LHXI,LGGFX, &
          data%ldx, data%lgrjac, data%lqgrad, data%lbreak, data%lp,     data%lxcp, data%lx0,  &
          data%lgx0, data%ldeltx, data%lbnd,   data%lwkstr, data%lsptrs, data%lselts, data%lindex,  &
          data%lswksp, data%lstagv, data%lstajc, data%liused, data%lfreec, data%lnnonz, data%lnonz2,  &
          data%lsymmd, data%lsymmh, data%lslgrp, data%lsvgrp, data%lgcolj, data%lvaljr, data%lsend,  &
          data%lnptrs, data%lnelts, data%lnndex, data%lnwksp, data%lnstgv, data%lnstjc, data%lniuse,  &
          data%lnfrec, data%lnnnon, data%lnnno2, data%lnsymd, data%lnsymh, data%lnlgrp, data%lnvgrp, &
          data%lngclj, data%lnvljr, data%lnqgrd, data%lnbrak, data%lnp,    data%lnbnd, &
          data%lnfxi,  data%lngxi,  data%lnguvl, data%lnhxi,  data%lnhuvl, data%lnggfx, &
          data%lndx, data%lngrjc, data%liwk2, data%lwk2, data%maxsin, data%ninvar, &
          data%ntype, data%nsets, data%maxsel, data%lstype, data%lsswtr, data%lssiwt, &
          data%lsiwtr, data%lswtra, data%lntype, data%lnswtr, data%lnsiwt, data%lniwtr, &
          data%lnwtra, data%lsiset, data%lssvse, data%lniset, data%lnsvse, RANGE, &
          data%IWORK(IWRK + 1),LIWORK,data%WRK,LWORK,IPRINT,IOUT,INFORM )
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

!  Initialize the performance counters and variables

      data%nc2of = 0
      data%nc2og = 0
      data%nc2oh = 0
      data%nc2cf = 0
      data%nc2cg = 0
      data%nc2ch = 0
      data%nhvpr = 0
      data%pnc = 0
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
 2000 FORMAT( /, ' ** Program USETUP: array length ', A6, ' too small.', &
              /, ' -- Miminimization abandoned.', &
              /, ' -- Increase the parameter ', A6, ' by at least ', I8, &
                 ' and restart.' )
 2010 FORMAT( /, ' ** Program USETUP: the problem includes general', &
                 ' constraints. Execution terminating ' )
 2020 FORMAT( /, ' ** Program USETUP: the problem uses no variables.', &
                 ' Execution terminating ' )
 2030 FORMAT( /, ' ** Program USETUP: the problem is vacuous.', &
                 ' Execution terminating ' )

!  End of USETUP.

      END
