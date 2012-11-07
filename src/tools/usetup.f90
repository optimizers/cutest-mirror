! THIS VERSION: CUTEST 1.0 - 04/11/2012 AT 12:30 GMT.

!-*-*-*-*-*-*-  C U T E S T    U S E T U P    S U B R O U T I N E  -*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 77 version originally released in CUTE, 30th October, 1991
!   fortran 2003 version released in CUTEst, 4th November 2012

      SUBROUTINE USETUP( data, input, iout, n, X, BL, BU, nmax )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: input, iout, n, nmax
      REAL ( KIND = wp ) :: X( nmax ), BL( nmax ), BU( nmax )

!  --------------------------------------------------------------
!  set up the input data for the unconstrained optimization tools
!  --------------------------------------------------------------

!  local variables

      INTEGER :: ialgor, iprint, inform, i, nslack
      LOGICAL :: fdgrad, debug
      REAL ( KIND = wp ), DIMENSION( 2 ) :: OBFBND
      CHARACTER ( LEN = 8 ) :: pname
      CHARACTER ( LEN = 10 ) :: chtemp
      EXTERNAL :: RANGE

      CALL CPU_TIME( data%sutime )
      data%iout2 = iout
      debug = .FALSE.
      debug = debug .AND. iout > 0
      iprint = 0
      IF ( debug ) iprint = 3

!  Input the problem dimensions.

      READ( input, 1001 ) data%n, data%ng, data%nel, data%ntotel, data%nvrels, &
                          data%nnza, data%ngpvlu, data%nepvlu, neltyp, ngrtyp
      n = data%n
      IF ( n <= 0 ) THEN
        CLOSE( input )
        IF ( iout > 0 ) WRITE( iout,                                           &
          "( /, ' ** Program USETUP: the problem uses no variables.',          &
         &       ' Execution terminating ' )" )
        STOP
      END IF
      IF ( data%ng <= 0 ) THEN
        CLOSE( input )
        IF ( iout > 0 ) WRITE( iout,                                           &
           "( /, ' ** Program USETUP: the problem is vacuous.',                &
          &      ' Execution terminating ' )" )
        STOP
      END IF
      IF ( n > nmax ) THEN
        CLOSE( input )
        IF ( iout > 0 ) THEN
          WRITE( iout, 2000 ) 'X', 'NMAX', n - nmax
          WRITE( iout, 2000 ) 'BL', 'NMAX', n - nmax
          WRITE( iout, 2000 ) 'BU', 'NMAX', n - nmax
        END IF
        STOP
      END IF

!  Input the problem type.

      READ( input, 1000 ) ialgor, pname

!  Set useful integer values.

      data%ng1 = data%ng + 1
      data%ngng = data%ng + data%ng
      data%nel1 = data%nel + 1

!  Partition the integer workspace.

      ISTADG = 0
      ISTGP = ISTADG + data%ng1
      ISTADA = ISTGP + data%ng1
      ISTAEV = ISTADA + data%ng1
      ISTEP = ISTAEV + data%nel1
      ITYPEG = ISTEP + data%nel1
      KNDOFC = ITYPEG + data%ng
      ITYPEE = KNDOFC + data%ng
      IELING = ITYPEE + data%nel
      IELVAR = IELING + data%ntotel
      ICNA = IELVAR + data%nvrels
      ISTADH = ICNA + data%nnza
      INTVAR = ISTADH + data%nel1
      IVAR = INTVAR + data%nel1
      ICALCF = IVAR + n
      ITYPEV = ICALCF + MAX( data%nel, data%ng )
      IWRK = ITYPEV + n
      data%liwork = liwk - IWRK

!  Ensure there is sufficient room.

      IF ( data%liwork < 0 ) THEN
        CLOSE( input )
        IF ( iout > 0 ) WRITE( iout, 2000 ) 'IWK', 'LIWK', - data%liwork
        STOP
      END IF

!  Partition the real workspace.

      A = 0
      B = A + data%nnza
      U = B + data%ng
      GPVALU = U + data%ng
      EPVALU = GPVALU + data%ngpvlu
      ESCALE = EPVALU + data%nepvlu
      GSCALE = ESCALE + data%ntotel
      VSCALE = GSCALE + data%ng
      GVALS = VSCALE + n
      XT = GVALS + 3 * data%ng
      DGRAD = XT + n
      Q = DGRAD + n
      FT = Q + n
      WRK = FT + data%ng
      data%lwork = lwk - WRK

!  Ensure there is sufficient room.

      IF ( data%lwork < 0 ) THEN
        CLOSE( input )
        IF ( iout > 0 ) WRITE( iout, 2000 ) 'WK', 'LWK', - data%lwork
        STOP
      END IF

!  Partition the logical workspace.

      intrep = 0
      gxeqx = intrep + data%nel
      data%lo = gxeqx + data%ngng

!  Ensure there is sufficient room.

      IF ( llogic < data%lo ) THEN
        CLOSE( input )
        IF ( iout > 0 ) WRITE( iout, 2000 ) 'LOGI', 'LLOGIC', data%lo - llogic
        STOP
      END IF

!  Partition the character workspace.

      GNAMES = 0
      VNAMES = GNAMES + data%ng
      data%ch = VNAMES + n

!  Ensure there is sufficient room.

      IF ( lchara < data%ch + 1 ) THEN
        CLOSE( input )
        IF ( iout > 0 ) WRITE( iout, 2000 ) 'CHA', 'LCHARA', data%ch + 1 -lchara
        STOP
      END IF

!  Record the lengths of arrays.

      data%lstadg = MAX( 1, data%ng1 )
      data%lstada = MAX( 1, data%ng1 )
      data%lstaev = MAX( 1, data%nel1 )
      data%lkndof = MAX( 1, data%ng )
      data%leling = MAX( 1, data%ntotel )
      data%lelvar = MAX( 1, data%nvrels )
      data%licna = MAX( 1, data%nnza )
      data%lstadh = MAX( 1, data%nel1 )
      data%lntvar = MAX( 1, data%nel1 )
      data%lcalcf = MAX( 1, data%nel, data%ng )
      data%lcalcg = MAX( 1, data%ng )
      data%la = MAX( 1, data%nnza )
      data%lb = MAX( 1, data%ng )
      data%lu = MAX( 1, data%ng )
      data%lescal = MAX( 1, data%ntotel )
      data%lgscal = MAX( 1, data%ng )
      data%lvscal = MAX( 1, n )
      data%lft = MAX( 1, data%ng )
      data%lgvals = MAX( 1, data%ng )
      data%lintre = MAX( 1, data%nel )
      data%lgxeqx = MAX( 1, data%ngng )
      data%lgpvlu = MAX( 1, data%ngpvlu )
      data%lepvlu = MAX( 1, data%nepvlu )
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

!  Print out problem data. input the number of variables, groups,
!  elements and the identity of the objective function group.

      IF ( ialgor == 2 ) THEN
        READ( input, 1002 ) nslack, data%nobjgr
      ELSE
        nslack = 0
      END IF
      IF ( debug ) WRITE( iout, 1100 ) pname, n, data%ng, data%nel
      data%pname = pname // '  '

!  Input the starting addresses of the elements in each group,
!  of the parameters used for each group and
!  of the nonzeros of the linear element in each group.

      READ( input, 1010 ) ( data%ISTADG( i ), i = 1, data%ng1 )
      IF ( debug ) WRITE( iout, 1110 ) 'ISTADG',                               &
        ( data%ISTADG( i ), i = 1, data%ng1 )
      READ( input, 1010 ) ( data%ISTGP( i ), i = 1, data%ng1 )
      IF ( debug ) WRITE( iout, 1110 ) 'ISTGP ',                               &
        ( data%ISTGP( i ), i = 1, data%ng1 )
      READ( input, 1010 ) ( data%ISTADA( i ), i = 1, data%ng1 )
      IF ( debug ) WRITE( iout, 1110 ) 'ISTADA',                               &
        ( data%ISTADA( i ), i = 1, data%ng1 )

!  Input the starting addresses of the variables and parameters
!  in each element.

      READ( input, 1010 ) ( data%ISTAEV( i ), i = 1, data%nel1 )
      IF ( debug ) WRITE( iout, 1110 ) 'ISTAEV',                               &
        ( data%ISTAEV( i ), i = 1, data%nel1 )
      READ( input, 1010 ) ( data%ISTEP( i ), i = 1, data%nel1 )
      IF ( debug ) WRITE( iout, 1110 ) 'ISTEP ',                               &
        ( data%ISTEP( i ), i = 1, data%nel1 )

!  Input the group type of each group

      READ( input, 1010 ) ( data%ITYPEG( i ), i = 1, data%ng )
      IF ( debug ) WRITE( iout, 1110 ) 'ITYPEG',                               &
        ( data%ITYPEG( i ), i = 1, data%ng )
      IF ( ialgor >= 2 ) THEN
        READ( input, 1010 )( data%KNDOFC( i ), i = 1, data%ng )
        IF ( debug ) WRITE( iout, 1110 ) 'KNDOFC',                             &
          ( data%KNDOFC( i ), i = 1, data%ng )
        DO 10 i = 1, data%ng
          IF ( ABS( data%KNDOFC( i ) ) >= 2 ) THEN
            CLOSE( input )
            IF ( iout > 0 ) WRITE( iout,                                       &
              "( /, ' ** Program USETUP: the problem includes general',        &
             &      ' constraints. Execution terminating ' )" )
            STOP
          END IF
   10   CONTINUE
      END IF

!  Input the element type of each element

      READ( input, 1010 ) ( data%ITYPEE( i ), i = 1, data%nel )
      IF ( debug ) WRITE( iout, 1110 ) 'ITYPEE',                               &
        ( data%ITYPEE( i ), i = 1, data%nel )

!  Input the number of internal variables for each element.

      READ( input, 1010 ) ( data%INTVAR( i ), i = 1, data%nel )
      IF ( debug ) WRITE( iout, 1110 ) 'INTVAR',                               &
        ( data%INTVAR( i ), i = 1, data%nel )

!  Input the identity of each individual element.

      READ( input, 1010 ) ( data%IELING( i ), i = 1, data%ntotel )
      IF ( debug ) WRITE( iout, 1110 ) 'IELING',                               &
        ( data%IELING( i ), i = 1, data%ntotel )

!  Input the variables in each group's elements.

      data%nvrels = data%ISTAEV( data%nel1 ) - 1
      READ( input, 1010 ) ( data%IELVAR( i ), i = 1, data%nvrels )
      IF ( debug ) WRITE( iout, 1110 ) 'IELVAR',                               &
        ( data%IELVAR( i ), i = 1, data%nvrels )

!  Input the column addresses of the nonzeros in each linear element.

      READ( input, 1010 ) ( data%ICNA( i ), i = 1, data%nnza )
      IF ( debug ) WRITE( iout, 1110 ) 'ICNA  ',                               &
        ( data%ICNA( i ), i = 1, data%nnza )

!  Input the values of the nonzeros in each linear element, the
!  constant term in each group, the lower and upper bounds on
!  the variables and the starting point for the minimization.

      READ( input, 1020 ) ( data%A( i ), i = 1, data%nnza )
      IF ( debug ) WRITE( iout, 1120 ) 'A     ',                               &
        ( data%A( i ), i = 1, data%nnza )
      READ( input, 1020 ) ( data%B( i ), i = 1, data%ng )
      IF ( debug ) WRITE( iout, 1120 ) 'B     ',                               &
        ( data%B( i ), i = 1, data%ng )
      IF ( ialgor <= 2 ) THEN
        READ( input, 1020 ) ( BL( i ), i = 1, n )
        IF ( debug ) WRITE( iout, 1120 ) 'BL    ', ( BL( i ), i = 1, n )
        READ( input, 1020 ) ( BU( i ), i = 1, n )
        IF ( debug ) WRITE( iout, 1120 ) 'BU    ', ( BU( i ), i = 1, n )
      ELSE

!  Use GVALS and FT as temporary storage for the constraint bounds.

        READ( input, 1020 ) ( BL( i ), i = 1, n ),                             &
          ( data%GVALS( i ), i = 1, data%ng )
        IF ( debug ) WRITE( iout, 1120 ) 'BL    ',                             &
          ( BL( i ), i = 1, n ), ( data%GVALS( i ), i = 1, data%ng )
        READ( input, 1020 ) ( BU( i ), i = 1, n ),                             &
          ( data%FT( i ), i = 1, data%ng )
        IF ( debug ) WRITE( iout, 1120 ) 'BU    ',                             &
          ( BU( i ), i = 1, n ), ( data%FT( i ), i = 1, data%ng )
      END IF
      READ( input, 1020 ) ( X( i ), i = 1, n )
      IF ( debug ) WRITE( iout, 1120 ) 'X     ', ( X( i ), i = 1, n )
      IF ( ialgor >= 2 ) THEN
        READ( input, 1020 )( data%U( i ), i = 1, data%ng )
        IF ( debug ) WRITE( iout, 1120 ) 'U     ',                             &
          ( data%U( i ), i = 1, data%ng )
      END IF

!  Input the parameters in each group.

      READ( input, 1020 ) ( data%GPVALU( i ), i = 1, data%ngpvlu )
      IF ( debug ) WRITE( iout, 1120 ) 'GPVALU',                               &
        ( data%GPVALU( i ), i = 1, data%ngpvlu )

!  Input the parameters in each individual element.

      READ( input, 1020 ) ( data%EPVALU( i ), i = 1, data%nepvlu )
      IF ( debug ) WRITE( iout, 1120 ) 'EPVALU',                               &
        ( data%EPVALU( i ), i = 1, data%nepvlu )

!  Input the scale factors for the nonlinear elements.

      READ( input, 1020 ) ( data%ESCALE( i ), i = 1, data%ntotel )
      IF ( debug ) WRITE( iout, 1120 ) 'ESCALE',                               &
        ( data%ESCALE( i ), i = 1, data%ntotel )

!  Input the scale factors for the groups.

      READ( input, 1020 ) ( data%GSCALE( i ), i = 1, data%ng )
      IF ( debug ) WRITE( iout, 1120 ) 'GSCALE',                               &
        ( data%GSCALE( i ), i = 1, data%ng )

!  Input the scale factors for the variables.

      READ( input, 1020 ) ( data%VSCALE( i ), i = 1, n )
      IF ( debug ) WRITE( iout, 1120 ) 'VSCALE',                               &
        ( data%VSCALE( i ), i = 1, n )

!  Input the lower and upper bounds on the objective function.

      READ( input, 1080 ) OBFBND( 1 ), OBFBND( 2 )
      IF ( debug ) WRITE( iout, 1180 ) 'OBFBND', OBFBND( 1 ), OBFBND( 2 )

!  Input a logical array which says whether an element has internal
!  varaiables.

      READ( input, 1030 ) ( data%INTREP( i ), i = 1, data%nel )
      IF ( debug ) WRITE( iout, 1130 ) 'INTREP',                               &
        ( data%INTREP( i ), i = 1, data%nel )

!  Input a logical array which says whether a group is trivial.

      READ( input, 1030 ) ( data%GXEQX( i ), i = 1, data%ng )
      IF ( debug ) WRITE( iout, 1130 ) 'GXEQX ',                               &
        ( data%GXEQX( i ), i = 1, data%ng )

!  Input the names given to the groups and to the variables.

      READ( input, 1040 ) ( data%GNAMES( i ), i = 1, data%ng )
      IF ( debug ) WRITE( iout, 1140 ) 'GNAMES',                               &
        ( data%GNAMES( i ), i = 1, data%ng )
      READ( input, 1040 ) ( data%VNAMES( i ), i = 1, n )
      IF ( debug ) WRITE( iout, 1140 ) 'VNAMES', ( data%VNAMES( i ), i = 1, n )

!  Dummy input for the names given to the element and group types.

      READ( input, 1040 ) ( chtemp, i = 1, neltyp )
      READ( input, 1040 ) ( chtemp, i = 1, ngrtyp )

!  Input the type of each variable.

      READ( input, 1010 ) ( data%ITYPEV( i ), i = 1, n )
      CLOSE( input )

      data%numvar = n

!  Partition the workspace arrays data%FUVALS, IWK and WK. Initialize
!  certain portions of IWK.

      data%firstg = .TRUE.
      fdgrad = .FALSE.
      CALL INITW( n, data%ng, data%nel, data%IELING, data%leling,           &
          data%ISTADG, data%lstadg, data%IELVAR, data%lelvar, data%ISTAEV,     &
          data%lstaev, data%INTVAR, data%lntvar, data%ISTADH, data%lstadh,     &
          data%ICNA, data%licna, data%ISTADA, data%lstada, data%ITYPEE,        &
          data%lintre, data%GXEQX, data%lgxeqx, data%INTREP, data%lintre,      &
          lfuval, data%altriv, .TRUE., fdgrad, data%lfxi, LGXI, LHXI, LGGFX,   &
          data%ldx, data%lgrjac, data%lqgrad, data%lbreak, data%lp,            &
          data%lxcp, data%lx0, data%lgx0, data%ldeltx, data%lbnd, data%lwkstr, &
          data%lsptrs, data%lselts, data%lindex, data%lswksp, data%lstagv,     &
          data%lstajc, data%liused, data%lfreec, data%lnnonz, data%lnonz2,     &
          data%lsymmd, data%lsymmh, data%lslgrp, data%lsvgrp, data%lgcolj,     &
          data%lvaljr, data%lsend, data%lnptrs, data%lnelts, data%lnndex,      &
          data%lnwksp, data%lnstgv, data%lnstjc, data%lniuse, data%lnfrec,     &
          data%lnnnon, data%lnnno2, data%lnsymd, data%lnsymh, data%lnlgrp,     &
          data%lnvgrp, data%lngclj, data%lnvljr, data%lnqgrd, data%lnbrak,     &
          data%lnp, data%lnbnd, data%lnfxi, data%lngxi, data%lnguvl,           &
          data%lnhxi, data%lnhuvl, data%lnggfx, data%lndx, data%lngrjc,        &
          data%liwk2, data%lwk2, data%maxsin, data%ninvar, data%ntype,         &
          data%nsets, data%maxsel, data%lstype, data%lsswtr, data%lssiwt,      &
          data%lsiwtr, data%lswtra, data%lntype, data%lnswtr, data%lnsiwt,     &
          data%lniwtr, data%lnwtra, data%lsiset, data%lssvse, data%lniset,     &
          data%lnsvse, RANGE, data%IWORK(IWRK + 1), liwork, data%WRK, lwork,   &
          iprint, iout, inform )
      IF ( inform /= 0 ) STOP

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

      CALL CPU_TIME( data%sttime )
      data%sutime = data%sttime - data%sutime
      RETURN

!  non-executable statements

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
 2000 FORMAT( /, ' ** Program USETUP: array length ', A, ' too small.',        &
              /, ' -- Miminimization abandoned.',                              &
              /, ' -- Increase the parameter ', A, ' by at least ', I0,        &
                 ' and restart.' )

!  End of subroutine USETUP

      END SUBROUTINE USETUP
