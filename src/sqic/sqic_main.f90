!     ( Last modified on 26 Feb 2013 at 13:30:00 )

  program sqic_main

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
!     Driver for running SQIC on CUTEst problems.
 
!     Derived from SQIC program sqicma.f90 written by 
!     Philip Gill and Elizabeth Wong
!     CUTEst evolution February 2013, Nick Gould
      
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! use snModulePrecision, only : ip
  use CUTEst_interface_double

  use ic21mxv,           only : crd2spr
  use sqModuleProb,      only : qpProb

  implicit none

  integer, parameter :: ip = kind( 1 ), rp = kind( 1.0D+0 )
  integer(ip), parameter :: iCutest = 55, iOut = 6, io_buffer = 11
  integer(ip), parameter :: iPrint = 9, iSumm  = 6, iSpecs = 4
  integer(ip)   :: INFO, n, m, nm, nnH, ncObj, neA, lenA, neH, lenH, nS, nInf
  integer(ip)   :: iObj, nName, status, alloc_stat
  real(rp)      :: Obj, ObjAdd, sInf
  real(rp)      :: CPU( 2 ), CALLS( 7 )
  character(8)  :: Names(1)
  character(10) :: Prob
  character(20) :: filename

  integer(ip),   pointer :: hs(:), hEtype(:), indA(:), locA(:)
  integer(ip),   pointer :: indH(:), locH(:), irow(:), jcol(:)
  real(rp),      pointer :: bl(:), bu(:), x(:), pi(:), rc(:), cObj(:)
  real(rp),      pointer :: valA(:), valH(:), zero(:), cval(:), b(:)
  logical,       pointer :: equation(:), linear(:)
  character(10), pointer :: vname(:), gname(:)
  type(qpProb)           :: QP

  open ( iCutest, file = 'OUTSDIF.d', form = 'FORMATTED', status = 'old' )
  rewind( iCutest )

  ! Get dimensions and allocate space.
  call CUTEST_cdimen( status, iCutest, n, m )
  if ( status /= 0 ) go to 910

  !-----------------------------------------------------------------------------
  ! Initial point and bounds ( n, m, x, bl, bu )
  ! Problem name ( Prob )
  ! Constraints

  if ( m > 0 ) then
     nm = n+m

     allocate ( bl(nm), bu(nm), cObj(n), x(nm), hs(nm), pi(m), rc(nm),         &
                hEtype(nm), b(m), zero(n), equation(m), linear(m),             &
                stat = alloc_stat )
     if ( alloc_stat /= 0 ) GO TO 990
     hEtype = 0
     zero   = 0.0
     x      = 0.0

     call CUTEST_csetup( status, iCutest, iOut, io_buffer,                     &
                         n, m, x(1:n), bl(1:n), bu(1:n),                       &
                         pi, bl(n+1:nm), bu(n+1:nm), equation, linear, 0, 0, 1 )
     if ( status /= 0 ) go to 910
     deallocate ( equation, linear )

     allocate ( vname(n), gname(m), stat = alloc_stat )
     if ( alloc_stat /= 0 ) GO TO 990

     call CUTEST_cnames( status, n, m, Prob, vname, gname )
     if ( status /= 0 ) go to 910
     deallocate ( vname, gname )

     call CUTEST_cdimsj( status, lenA )
     if ( status /= 0 ) go to 910
     lenA = lenA - n
     allocate ( cval(lenA), irow(lenA), jcol(lenA), stat = alloc_stat )
     if ( alloc_stat /= 0 ) GO TO 990

     call CUTEST_ccfsg( status, n, m, zero, b, neA, lenA, cval, jcol, irow,    &
                        .true. )
     if ( status /= 0 ) go to 910

     bl(n+1:nm) = bl(n+1:nm) - b
     bu(n+1:nm) = bu(n+1:nm) - b
     deallocate ( b )

     ! Convert coordinate form to sparse-by-col form.
     allocate ( indA(neA), locA(n+1), valA(neA), stat = alloc_stat )
     if ( alloc_stat /= 0 ) GO TO 990
     call crd2spr ( n, neA, cval, irow, jcol, valA, indA, locA )
     deallocate ( irow, jcol, cval )

     ! Objective
     ! cofg returns ObjAdd and cObj.
     call CUTEST_cofg( status, n, zero, ObjAdd, cObj, .true. )
     if ( status /= 0 ) go to 910

     ! Hessian of the objective
     ! Stored in the module variables neH, valH, locH, indH.
     call CUTEST_cdimsh ( status, lenH )
     if ( status /= 0 ) go to 910

     if ( lenH > 0 ) then
        allocate ( cval(lenH), irow(lenH), jcol(lenH), stat = alloc_stat )
        if ( alloc_stat /= 0 ) GO TO 990
        call CUTEST_cish ( status, n, zero, 0, neH, lenH, cval, irow, jcol )
        if ( status /= 0 ) go to 910

        deallocate ( zero )

        nnH = maxval ( jcol(1:neH) )

        allocate ( indH(neH), locH(n+1), valH(neH), stat = alloc_stat )
        if ( alloc_stat /= 0 ) GO TO 990
        call crd2spr ( n, neH, cval, irow, jcol, valH, indH, locH )
        deallocate ( irow, jcol, cval )
     else
        nnH = 0
     end if

  else

     m  = 1
     nm = n+m

     allocate ( bl(nm), bu(nm), cObj(n), x(nm), hs(nm), pi(m), rc(nm),         &
                hEtype(nm), b(m), zero(n), stat = alloc_stat )
     if ( alloc_stat /= 0 ) GO TO 990
     hEtype = 0
     zero   = 0.0
     x      = 0.0

     call CUTEST_usetup ( status, iCutest, iOut, io_buffer,                    &
                          n, x(1:n), bl(1:n), bu(1:n) )
     if ( status /= 0 ) go to 910

     allocate ( vname(n), stat = alloc_stat )
     if ( alloc_stat /= 0 ) GO TO 990
     call CUTEST_unames ( status, n, Prob, vname )
     if ( status /= 0 ) go to 910
     deallocate ( vname )

     ! dummy row in A
     neA = 1
     allocate ( indA(neA), locA(n+1), valA(neA), stat = alloc_stat )
     if ( alloc_stat /= 0 ) GO TO 990
     indA(1)   = 1
     valA(1)   = 1.0_rp
     locA      = 1
     locA(n+1) = 2

     x(n+1)    =  0.0
     bl(n+1)   = -1.0d+20
     bu(n+1)   =  1.0d+20

     ! Objective
     ! uofg returns ObjAdd and cObj.
     call CUTEST_uofg( status, n, zero, ObjAdd, cObj, .true. )
     if ( status /= 0 ) go to 910

     ! Hessian of the objective
     call CUTEST_udimsh ( status, lenH )
     if ( status /= 0 ) go to 910

     if ( lenH > 0 ) then
        allocate ( cval(lenH), irow(lenH), jcol(lenH), stat = alloc_stat )
        if ( alloc_stat /= 0 ) GO TO 990
        call CUTEST_ush( status, n, zero, neH, lenH, cval, irow, jcol )
        if ( status /= 0 ) go to 910

        deallocate ( zero )

        nnH = maxval ( jcol(1:neH) )

        allocate ( indH(neH), locH(n+1), valH(neH), stat = alloc_stat )
        if ( alloc_stat /= 0 ) GO TO 990
        call crd2spr ( n, neH, cval, irow, jcol, valH, indH, locH )
        deallocate ( irow, jcol, cval )

     else
        nnH = 0
     end if
  end if

! write(6,*) ' m, n ', m, n, Prob

  !-----------------------------------------------------------------------------
  ! Ok, we're done with CUTEst stuff.
  !-----------------------------------------------------------------------------

  filename = trim(Prob)//'.out'
  open ( iPrint, file=filename, status='unknown' )
  open ( iSpecs, file='SQIC.SPC', form='formatted', status='old' )

  iObj  = 0
  ncObj = n
  nName = 1
  hs    = 0

  call qp%begin ( iPrint, iSumm )
  call qp%specs ( iSpecs, INFO )

  call qp%load ( Prob, m, n, nnH, iObj, ObjAdd,                                &
                 neA, indA, locA, valA, bl, bu, ncObj, cObj,                   &
                 nName, Names, hEtype, hs, x, pi, rc, neH, indH, locH, valH )

  call qp%solve ( 'Cold', INFO, nS, nInf, sInf, Obj )
!  call qp%solve ( 'Warm', INFO, nS, nInf, sInf, Obj )
!  call qp%solve ( 'Hot', INFO, nS, nInf, sInf, Obj )
  call qp%end

  call CUTEST_creport( status, CALLS, CPU )
  WRITE ( iOut, 2000 ) Prob, n, m, CALLS( 1 ), CALLS( 2 ),                     &
                      CALLS( 5 ), CALLS( 6 ), info, Obj, CPU( 1 ), CPU( 2 )

  deallocate ( bl, bu, cObj, x, pi, rc )
  deallocate ( hEtype, hs )
  deallocate ( indA, locA, valA )

  if ( nnH > 0 ) deallocate ( indH, locH, valH )

  close ( iSpecs )
  close ( iPrint )
  close ( iCutest )

  call CUTEST_cterminate( status )
  stop

  910 CONTINUE
  WRITE( iOut, "( ' CUTEst error, status = ', i0, ', stopping' )") status
  stop

  990 CONTINUE
  WRITE( iOut, "( ' Allocation error, status = ', i0 )" ) status
  stop

 2000 FORMAT( /, 24('*'), ' CUTEst statistics ', 24('*') //                    &
          ,' Package used            :  SQIC',    /                            &
          ,' Problem                 :  ', A10,    /                           &
          ,' # variables             =      ', I10 /                           &
          ,' # constraints           =      ', I10 /                           &
          ,' # objective functions   =        ', F8.2 /                        &
          ,' # objective gradients   =        ', F8.2 /                        &
          ,' # constraints functions =        ', F8.2 /                        &
          ,' # constraints gradients =        ', F8.2 /                        &
          ,' Exit code               =      ', I10 /                           &
          ,' Final f                 = ', E15.7 /                              &
          ,' Set up time             =      ', 0P, F10.2, ' seconds' /         &
          ,' Solve time              =      ', 0P, F10.2, ' seconds' //        &
           66('*') / )

  end program sqic_main
