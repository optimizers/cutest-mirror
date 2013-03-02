!     ( Last modified on 2 Mar 2013 at 13:30:00 )

!  Dummy SQIC modules for testing sqic_main interface to CUTEst
!  Nick Gould, 2nd March 2013

module ic21mxv
 public  :: crd2spr
  integer, parameter :: ip = kind( 1 ), rp = kind( 1.0D+0 )
contains
 subroutine crd2spr ( n, ne, cval, irow, jcol, val, ind, loc )
  implicit none

  integer(ip), parameter :: ip = kind( 1 ), rp = kind( 1.0D+0 )
  integer(ip), intent(in)  :: n, ne, irow(ne), jcol(ne)
  real(rp),    intent(in)  :: cval(ne)
  integer(ip), intent(out) :: ind(ne), loc(n+1)
  real(rp),    intent(out) :: val(ne)
  !=============================================================================
  ! Converts coordinate form (cval, irow, jcol) into sparse-by-column format in
  ! val, ind, loc.
  !
  ! Philip Gill and Elizabeth Wong
  ! 15 Feb 2010: First version.
  ! ic21mxv.f90 from sqic distribution.
  !=============================================================================
  integer(ip) :: i, j, k, ii

  ! First count up elements in each column and store in loc.
  loc      = 0
  do k = 1, ne
     j = jcol(k)
     i = irow(k)
     loc(j) = loc(j) + 1
  end do
  loc(n+1) = ne + 1

  ! Set the column pointers.
  do k = n, 1, -1
     loc(k) = loc(k+1) - loc(k)
  end do

  ! Put the row indices and values in and let loc(j) track the position of the
  ! elements for the jth column.
  ind = 0
  val = 0.0
  do k = 1, ne
     j  = jcol(k)
     i  = irow(k)
     ii = loc(j)

     val(ii) = cval(k)
     ind(ii) = i
     loc(j)  = ii + 1
  end do

  ! Reset the column pointers.
  loc(2:n) = loc(1:n-1)
  loc(1)   = 1

 end subroutine crd2spr
end module ic21mxv

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! slightly modified header from snWork.f90
! Philip Gill and Elizabeth Wong
! 15 Feb 2010: First version.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module snModuleWork
  implicit none
  private
  public  :: snWork, chkStart
  private :: dfltCmn, qpDflt, snDflt, allocCmn, qpAlloc, snAlloc, allocInt,    &
             allocReal, deallocInt, deallocReal, deallocate, undefine, copy

  integer, parameter :: ip = kind( 1 ), rp = kind( 1.0D+0 )
  integer(ip),  parameter :: blkQP = 10, varQP = 11, regQP = 12, dualQP = 13,  &
                             pdrQP = 14
  type snWork
     real(rp)    :: tolFP, tolQP, tolNLP, tolCG, tolx, tolCon, tolpiv, tolrow, &
                    tCrash, Utol1m, Utol2m, tolSwp, tolFac, tolUpd, LUspace,   &
                    infBnd, bigFx, bigdx, epsrf, fdint(2), xdlim, vilim,       &
                    etarg, wolfeG, Hcndbd, Zcndbd, wtInf0, xPen0, wtMax,       &
                    scltol, Aijtol, bStruc(2)
     integer(ip) :: lvlDer, lvlDerA, lvlInf, lvlHes, lvlPiv, lvlPPm, lvlPre,   &
                    lvlScl, lvlSch, lvlSrt, lvlSys, lvlTim, lvlVer
     integer(ip) :: maxR, maxS, mQNmod, QPslvr, Emode, kchk, kfac,             &
                    ksav, klog, kSumm, kDegen, kReset, lprPrm, lprSch,         &
                    lprScl, lprSol, lprDbg, minmax, mFlush, mSkip, iCrash,     &
                    itnlim, mMajor, mMinor, MjrPrt, MnrPrt, nParPr, mNewSB,    &
                    cgItmx, ObjRow, DenJ, MaxErr, MaxList, nProb, verify(6),   &
                    stkyOp, qpStat, npStat, iBack, iDump, iLoadB, iMPS,        &
                    iNewB, iInsrt, iOldB, iPnch, iReprt, iSoln, maxm, maxn,    &
                    maxne, DerOpt
     character(8) :: mProb, mObj, mRhs, mRng, mBnd
     integer(ip) :: LUparm(30)
     real(rp) :: parmLU(30)
     integer(ip) :: iStdi, iStdo, iSpecs, iPrint, iSumm, iPrinx, iSummx
     integer(ip) :: nnObj, nnCon, nnJac, lenR, negCon, ngObj, ngQP, nkx,       &
                    lenLU, lenMemI, lenMemR
     integer(ip) :: lvlDif
     integer(ip) :: iError, iALONE, iPage1, iPage2, PrintP, PrintS, linesP,    &
                    linesS, MnrHdP, MjrHdP, MnrHdS, MjrHdS, minimz, kObj,      &
                    nFac, nBfac, LUreq, LUitn, LUmod, QPmode, PreCon, gotFD,   &
                    gotG, gotGl, nfCon(4), nfObj(4), maxvi, HQNType, gotH,     &
                    nHess, nQPHx, nNPHx, eigH, gotFac, gotHes, gotScl, nName,  &
                    nQNmod, nSkip, nFlush, cgItn, cgItns
     real(rp) :: eps, eps0, eps1, eps2, eps3, eps4, eps5, tolDpp, tolDcp,      &
                 tolDup, toldj(3), tolDrp, tolDdp, Umin, condZ, yObj
     integer(ip) :: itn, nMajor, iExit, jbInf1, jdInf1, jbInf, jdInf
     real(rp)    :: ObjTru, piNorm, xNorm, wtInf, Binf1, Dinf1, Binf, Dinf,    &
                    sInf, fLin, fObj, PenNrm
     integer(ip) :: numt(10)
     real(rp)    :: tsum(10), tlast(10)
     integer(ip) :: itnfix, nfix(2)
     real(rp) :: featol, tolx0, tolinc
     character, pointer :: Names(:) * 8
     integer(ip), pointer :: kBS(:), keynm(:), kx(:), kxN(:), locJ(:),         &
                             indJ(:), locG(:), indG(:), nGlin(:), iGfun(:),    &
                             jGvar(:), hs(:), hElast(:), hfeas(:), hEstate(:)
     real(rp), pointer :: Ascale(:), bl(:), bu(:), blwrk(:), buwrk(:),         &
                          blBS(:),  buBS(:), dx(:), F(:), fCon(:), fCon1(:),   &
                          fCon2(:), Fmul(:), Fv(:), Fx(:), G(:), gBS(:),       &
                          gCon(:), gCon1(:), gCon2(:), gConu(:), gd(:),        &
                          gObj(:), gObj1(:), gObj2(:), gObju(:), gQP(:),       &
                          Gsave(:), Hdx(:), Hd(:), Ajac(:), pBS(:), pi(:),     &
                          QPrhs(:), R(:), rc(:), rg(:), rg2(:), Ux(:),         &
                          U0(:), x0(:), x(:), x1(:), xBS(:), xN(:),            &
                          xPen(:), xQP(:), xQP0(:), xscale(:), yCon(:),        &
                          yCon1(:), yCon2(:), dyCon(:)
     integer(ip),  pointer :: iy(:), iy1(:), iy2(:)
     real(rp),     pointer :: y(:), y1(:), y2(:)
     real(rp), pointer :: t1(:), r1(:), r2(:), s1(:), s2(:), s3(:)
     integer(ip),  pointer :: p(:), q(:), lenc(:), lenri(:), locc(:), locr(:), &
                              iploc(:), iqloc(:), indc(:), indr(:)
     real(rp),    pointer :: aLU(:)
     integer(ip)  :: itrLim, schLim, dbgLvl, lvlInd, itnChk, resChk, kktRef,   &
                     sclItn, iLU, dbgErr, vertex, rtry, nkFac, gotBlk
     integer(ip)  :: iAdd, iBnd, iOpt, iRmv, iSwp
     real(rp)     :: mu, rho, blkBnd, refBnd, schBnd, Anrm, Hnrm, rsdtol,      &
                     timeLim, tolInd, toluv, zero
     integer(ip), pointer :: hrtype(:), ztype(:)
     real(rp),    pointer :: pB(:), qB(:), uB(:), vB(:), Hvs(:), Avs(:), rsd(:)

     ! For external LU solvers:
     integer(ip) :: LUpiv
     real(rp)    :: lutol

   contains

     procedure, public,  pass :: init
     procedure, public,  pass :: undefine
     procedure, public,  pass :: copy
     procedure, public,  pass :: chkStart
     procedure, public,  pass :: qpAlloc
     procedure, public,  pass :: snAlloc
     procedure, private, pass :: allocCmn
     procedure, public,  pass :: deallocate
     procedure, public,  pass :: qpDflt
     procedure, public,  pass :: snDflt
     procedure, private, pass :: dfltCmn

  end type snWork

contains

  subroutine chkStart ( w, iExit, Start, lvlSrt, Errors )

    class(snWork), target        :: w
    character*(*)                :: Start
    integer(ip),   intent(inout) :: lvlSrt
    integer(ip),   intent(out)   :: iExit, Errors
  end subroutine chkStart

  subroutine init ( w )
    class(snWork) :: w
  end subroutine init

  subroutine dfltCmn ( w, m, n, linear, QP, nonlin )
    class(snWork), target     :: w
    logical,       intent(in) :: linear, QP, nonlin
    integer(ip),   intent(in) :: m, n

  end subroutine dfltCmn

  subroutine qpDflt ( w, m, n, nnH )
    class(snWork), target     :: w
    integer(ip),   intent(in) :: m, n, nnH
  end subroutine qpDflt

  subroutine snDflt ( w, m, n, ne, nnH )
    class(snWork), target     :: w
    integer(ip),   intent(in) :: m, n, ne, nnH
  end subroutine snDflt

  subroutine allocCmn ( w, m, n, nb, nnH, maxS, mBS, Err )
    class(snWork), target      :: w
    integer(ip),   intent(in)  :: m, n, nb, nnH, maxS, mBS
    integer(ip),   intent(out) :: Err
  end subroutine allocCmn

  subroutine qpAlloc ( w, iExit, m, n, ne, nnH )
    class(snWork), target      :: w
    integer(ip),   intent(in)  :: m, n, ne, nnH
    integer(ip),   intent(out) :: iExit
  end subroutine qpAlloc

  subroutine snAlloc ( w, iExit, m, n, ne, nnH, nF, neG )
    class(snWork), target      :: w
    integer(ip),   intent(in)  :: m, n, ne, nnH, nF, neG
    integer(ip),   intent(out) :: iExit
  end subroutine snAlloc

  subroutine deallocate ( w )
    class(snWork) :: w
  end subroutine deallocate

  subroutine undefine ( w )
    class(snWork) :: w
  end subroutine undefine

  subroutine copy ( w, w1 )
    class(snWork) :: w, w1
  end subroutine copy

  subroutine allocInt( Error, Name, iArray, lenReq, w )
    character,     intent(in)    :: Name*(*)
    integer(ip),   intent(in)    :: lenReq
    integer,       intent(inout) :: Error
    integer(ip),   pointer       :: iArray(:)
    class(snWork), target        :: w
  end subroutine allocInt

  subroutine allocReal( Error, Name, rArray, lenReq, w )
    integer,       intent(inout) :: Error
    character,     intent(in)    :: Name*(*)
    integer(ip),   intent(in)    :: lenReq
    real(rp),      pointer       :: rArray(:)
    class(snWork), target        :: w
  end subroutine allocReal

  subroutine deallocInt ( Error, Name, iArray, w )
    integer,       intent(inout) :: Error
    character,     intent(in)    :: Name*(*)
    integer(ip),   pointer       :: iArray(:)
    class(snWork), target        :: w
  end subroutine deallocInt

  subroutine deallocReal  ( Error, Name, rArray, w )
    integer,       intent(inout) :: Error
    character,     intent(in)    :: Name*(*)
    real(rp),      pointer       :: rArray(:)
    class(snWork), target        :: w
  end subroutine deallocReal

end module snModuleWork

!=============================================================================
! slightly modified header from ic45blk.f90
! Philip Gill and Elizabeth Wong
! 20 Dec 2011: First version with type extensions.
!=============================================================================

module ic45blk
  use  snModuleWork,      only : snWork
  implicit none
  private
  public  :: blkMx
  private :: q4fac, q4chk, q4sol, q4add, q4del, q4swp, q4trash, q4aug, q4brsz
  private :: sprMx, q2spr, q2add, q2del, q2trsh, q2Aprd

  integer, parameter :: ip = kind( 1 ), rp = kind( 1.0D+0 )

  type sprMx
     integer(ip), pointer :: ind(:) => null(), loc(:) => null(),               &
                             ptr(:) => null()
     real(rp),    pointer :: val(:) => null()

     integer(ip) :: se = 0, nc = 0, sc = 0

   contains
     procedure :: init  => q2spr
     procedure :: add   => q2add
     procedure :: del   => q2del
     procedure :: trash => q2trsh
     procedure :: q2Aprd
  end type sprMx

  type blkmx
     type(sprMx) :: Y, Z, V
     logical     :: done = .false.
     integer(ip) :: nB0  = 0, nB1  = 0, nK   = 0, nblk = 0, nV   = 0,          &
                    mxsc = 0, npos = 0, nneg = 0, rank = 0
     real(rp) :: condB = 0.0


     integer(ip), pointer :: sB0(:)  => null(), sB1(:)  => null(),             &
                             ind0(:) => null(), ind1(:) => null()
     real(rp),    pointer :: rhs(:) => null()
     real(rp), pointer :: wy(:,:) => null(), wz(:,:) => null(), &
                          wd(:,:) => null(), wx(:,:) => null()
   contains
     procedure :: fac   => q4fac
     procedure :: add   => q4add
     procedure :: del   => q4del
     procedure :: swp   => q4swp
     procedure :: sol   => q4sol
     procedure :: chk   => q4chk
     procedure :: trash => q4trash

  end type blkMx

contains

  subroutine q4fac ( blkB )
  class(blkMx), target :: blkB
  end subroutine q4fac

  subroutine q4chk ( blkB )
  class(blkMx), target :: blkB
  end subroutine q4chk

  subroutine q4sol ( blkB )
  class(blkMx), target :: blkB
  end subroutine q4sol

  subroutine q4add ( blkB )
  class(blkMx), target :: blkB
  end subroutine q4add

  subroutine q4del ( blkB )
  class(blkMx), target :: blkB
  end subroutine q4del

  subroutine q4swp ( blkB )
  class(blkMx), target :: blkB
  end subroutine q4swp

  subroutine q4aug ( blkB )
  class(blkMx), target :: blkB
  end subroutine q4aug

  subroutine q4brsz ( blkB )
  class(blkMx), target :: blkB
  end subroutine q4brsz

  subroutine q4trash ( blkB )
  class(blkMx), target :: blkB
  end subroutine q4trash

  subroutine q2spr ( Amtx )
  class(sprMx) :: Amtx
  end subroutine q2spr

  subroutine q2add ( Amtx )
  class(sprMx) :: Amtx
  end subroutine q2add

  subroutine q2del ( Amtx )
  class(sprMx) :: Amtx
  end subroutine q2del

  subroutine q2trsh ( Amtx )
  class(sprMx) :: Amtx
  end subroutine q2trsh

  subroutine q2Aprd ( Amtx )
  class(sprMx) :: Amtx
  end subroutine q2Aprd

end module ic45blk

!=============================================================================
! slightly modified header from snInterfaces.f90
! Philip Gill and Elizabeth Wong
! 20 Dec 2011: First version with type extensions.
!=============================================================================

module snInterfaces
  implicit none
  public
  integer, parameter :: ip = kind( 1 ), rp = kind( 1.0D+0 )

  interface

     subroutine iusrfunA ( Status, n, x, needF, nF, F, needG, lenG, G )
       integer, parameter :: ip = kind( 1 ), rp = kind( 1.0D+0 )
       integer(ip)  :: Status, n, needF, nF, needG, lenG
       real(rp)     :: F(nF), G(lenG), x(n)
     end subroutine iusrfunA

     subroutine iusrfunC ( mode, nnObj, nnCon, nnJac, nnL, neJac, x,           &
                           fObj, gObj, fCon, gCon, nState )
       integer, parameter :: ip = kind( 1 ), rp = kind( 1.0D+0 )
       integer(ip), intent(in)    :: nnObj, nnCon, nnJac, nnL, neJac, nState
       integer(ip), intent(inout) :: mode

       real(rp), intent(out) :: x(nnL), fObj, gObj(nnObj),                     &
                                fCon(nnCon), gCon(nnObj)
     end subroutine iusrfunC

     subroutine ifuncon ( mode, nnCon, nnJac, neJac, x, fCon, gCon, nState )
       integer, parameter :: ip = kind( 1 ), rp = kind( 1.0D+0 )
       integer(ip) :: mode, nnCon, nnJac, neJac, nState
       real(rp)    :: x(nnJac), fCon(nnCon), gCon(neJac)
     end subroutine ifuncon

     subroutine ifunobj ( mode, nnObj, x, fObj, gObj, nState )
       integer, parameter :: ip = kind( 1 ), rp = kind( 1.0D+0 )
       integer(ip) :: mode, nnObj, nState
       real(rp)    :: x(nnObj), fObj, gObj(nnObj)
     end subroutine ifunobj

     subroutine iusrHx ( nnH, x, Hx, jcol, State )
       integer, parameter :: ip = kind( 1 ), rp = kind( 1.0D+0 )
       integer(ip), intent(in)  :: nnH, jcol, State
       real(rp),    intent(in)  :: x(nnH)
       real(rp),    intent(out) :: Hx(nnH)
     end subroutine iusrHx

  end interface

end module snInterfaces

!=============================================================================
! slightly modified header from snInfo.f90
! Philip Gill and Elizabeth Wong
! 20 Dec 2011: First version with type extensions.
!=============================================================================

module snModuleInfo
! use  snModulePrecision, only : ip, rp
  use  snModuleWork,      only : snWork
! use  snModuleIO,        only : snPRNT
! use  snModuleMatrix,    only : snAj, snAij, snH
  use  snInterfaces

  implicit none

  public  :: snInfo, snInfoA, snInfoC, sqInfo

! integer, parameter :: ip = kind( 1 ), rp = kind( 1.0D+0 )

  type, abstract :: snInfo
     character(8) :: ProbName
     integer(ip)  :: n
     real(rp)     :: ObjAdd
   contains
     procedure(check_int), public, deferred, pass  :: check
  end type snInfo

  interface
     subroutine check_int ( info )
       use snModuleWork,      only : snWork
       import snInfo
       class(snInfo), target        :: info
     end subroutine check_int
  end interface

  type, extends(snInfo) :: snInfoA
     integer(ip)  :: nF, nxname, nFname, ObjRow
     character(8), pointer :: xnames(:), Fnames(:)
     integer(ip),  pointer :: xstate(:), Fstate(:)
     real(rp),     pointer :: xlow(:), xupp(:), Flow(:), Fupp(:),              &
                              x(:), F(:), xmul(:), Fmul(:)
     procedure(iusrfunA), nopass, pointer :: userfun => null()
   contains
     procedure, public, pass  :: load  => loadA
     procedure, public, pass  :: check => checkA
  end type snInfoA

  type, extends(snInfo) :: snInfoC
     integer(ip)  :: m, ne, nNames, nnCon, nnObj, nnJac, iObj
     character(8), pointer :: Names(:)
     integer(ip),  pointer :: hs(:)
     real(rp),     pointer :: bl(:), bu(:), x(:), pi(:), rc(:)
     procedure(iusrfunC), nopass, pointer :: userfun => null()
     procedure(ifuncon),  nopass, pointer :: funcon => null()
     procedure(ifunobj),  nopass, pointer :: funobj => null()
   contains
     procedure, public, pass  :: load  => loadC
     procedure, public, pass  :: check => checkC
  end type snInfoC

  type, extends(snInfo) :: sqInfo
     integer(ip)           :: m, ncObj, iObj, nNames
     character(8), pointer :: Names(:)
     integer(ip),  pointer :: hEtype(:), hs(:)
     real(rp),     pointer :: bl(:), bu(:), cObj(:), x(:), pi(:), rc(:)
   contains
     procedure, public, pass   :: load  => loadQ
     procedure, public, pass   :: loadH => loadQH
     procedure, public, pass   :: check => checkQ
  end type sqInfo

contains

  subroutine loadA ( info )
    class(snInfoA), target     :: info
  end subroutine loadA

  subroutine loadC (info  )
    class(snInfoC), target     :: info
  end subroutine loadC

  subroutine loadQ ( info )
    class(sqInfo), target     :: info
  end subroutine loadQ

  subroutine loadQH ( info )
    class(sqInfo), target     :: info
  end subroutine loadQH

  subroutine checkA ( info )
   class(snInfoA), target     :: info
  end subroutine checkA

  subroutine checkC ( info )
    class(snInfoC), target     :: info
  end subroutine checkC

  subroutine checkQ ( info )
    class(sqInfo), target     :: info
  end subroutine checkQ

end module snModuleInfo

!=============================================================================
! slightly modified header from snProb.f90
! Philip Gill and Elizabeth Wong
! 20 Dec 2011: First version with type extensions.
!=============================================================================

module snModuleProbs
  use  snModuleWork,      only : snWork
  use  ic45blk,           only : blkMx

  implicit none

  public  :: snProbs
  private :: getTitle, begin, end, specs, &
             set, seti, setr, get, getc, geti, getr

  integer, parameter :: ip = kind( 1 ), rp = kind( 1.0D+0 )

  type, abstract :: snProbs
     type(snWork) :: wrk
     type(blkMx)  :: blkB
   contains
     procedure, private, nopass :: getTitle
     procedure, public,  pass   :: begin
     procedure, public,  pass   :: end
     procedure, public,  pass   :: specs
     procedure, public,  pass   :: set
     procedure, public,  pass   :: seti
     procedure, public,  pass   :: setr
     procedure, public,  pass   :: get
     procedure, public,  pass   :: getc
     procedure, public,  pass   :: geti
     procedure, public,  pass   :: getr
     procedure(solve_int), public,  deferred, pass :: solve
     procedure(print_int), public,  deferred, pass :: print
     procedure(dflts_int), private, deferred, pass :: defaults
     procedure(alloc_int), private, deferred, pass :: alloc
  end type snProbs

  interface
     subroutine solve_int ( prob, Start, INFO, nS, nInf, sInf, Obj )
       import snProbs
       integer, parameter :: ip = kind( 1 ), rp = kind( 1.0D+0 )
       class(snProbs)               :: prob
       character*(*), intent(in)    :: Start
       integer(ip),   intent(inout) :: nInf
       integer(ip),   intent(out)   :: INFO, nS
       real(rp),      intent(out)   :: sInf, Obj
     end subroutine solve_int

     subroutine print_int ( prob )
       import snProbs
       class(snProbs), target :: prob
     end subroutine print_int

     subroutine dflts_int ( prob )
       import snProbs
       class(snProbs) :: prob
     end subroutine dflts_int

     subroutine alloc_int ( prob, inform )
       import snProbs
       integer, parameter :: ip = kind( 1 ), rp = kind( 1.0D+0 )
       class(snProbs)           :: prob
       integer(ip), intent(out) :: inform
     end subroutine alloc_int

  end interface

contains

  subroutine getTitle( title )
    character(30), intent(out) :: title
  end subroutine getTitle

  subroutine begin ( prob, iPrint, iSumm )
    class(snProbs), target     :: prob
    integer(ip),    intent(in) :: iPrint, iSumm
  end subroutine begin

  subroutine end ( prob )
    class(snProbs) :: prob
  end subroutine end

  subroutine specs ( prob, iSpecs, iExit )
    class(snProbs), target      :: prob
    integer(ip),    intent(in)  :: iSpecs
    integer,        intent(out) :: iExit
  end subroutine specs

  subroutine set ( prob, buffer, iPrint, iSumm, Errors )
    class(snProbs), target        :: prob
    character*(*),  intent(in)    :: buffer
    integer(ip),    intent(in)    :: iPrint, iSumm
    integer(ip),    intent(inout) :: Errors
  end subroutine set

  subroutine seti( prob, buffer, ivalue, iPrint, iSumm, Errors )
    class(snProbs), target        :: prob
    character*(*),  intent(in)    :: buffer
    integer(ip),    intent(in)    :: iPrint, iSumm
    integer(ip),    intent(inout) :: Errors
    integer(ip),    intent(in)    :: ivalue
  end subroutine Seti

  subroutine setr ( prob, buffer, rvalue, iPrint, iSumm, Errors )
    class(snProbs), target        :: prob
    character*(*),  intent(in)    :: buffer
    integer(ip),    intent(in)    :: iPrint, iSumm
    integer(ip),    intent(inout) :: Errors
    real(rp),       intent(in)    :: rvalue
  end subroutine setr

  function get  ( prob, buffer, Errors ) result(iget)
    class(snProbs), target        :: prob
    character*(*),  intent(in)    :: buffer
    integer(ip),    intent(inout) :: Errors
    integer(ip)                   :: iget
  end function get

  subroutine getc( prob, buffer, cvalue, Errors )
    class(snProbs), target        :: prob
    character,      intent(in)    :: buffer*(*)
    character(8),   intent(out)   :: cvalue
    integer(ip),    intent(inout) :: Errors
  end subroutine getc

  subroutine geti( prob, buffer, ivalue, Errors )
    class(snProbs), target        :: prob
    character*(*),  intent(in)    :: buffer
    integer(ip),    intent(out)   :: ivalue
    integer(ip),    intent(inout) :: Errors
  end subroutine geti

  subroutine getr( prob, buffer, rvalue, Errors )
    class(snProbs), target        :: prob
    character*(*),  intent(in)    :: buffer
    real(rp),       intent(out)   :: rvalue
    integer(ip),    intent(inout) :: Errors
  end subroutine getr

end module snModuleProbs

!=============================================================================
! slightly modified header from qpProb.f90
! Philip Gill and Elizabeth Wong
! 20 Dec 2011: First version with type extensions.
!=============================================================================

module sqModuleProb
  use  snInterfaces,      only : iusrHx
  use  snModuleProbs,     only : snProbs
  use  snModuleInfo,      only : sqInfo

  implicit none

  public  :: qpProb
  private :: sqTitle, loadQ, loadQH, printQ

  integer, parameter :: ip = kind( 1 ), rp = kind( 1.0D+0 )

  type, extends(snProbs) :: qpProb
     type(sqInfo) :: info
   contains
     procedure, private, nopass :: getTitle => sqTitle
     procedure, private, pass   :: loadQ
     procedure, private, pass   :: loadQH
     generic                    :: load => loadQ, loadQH
     procedure, public,  pass   :: solve
     procedure, public,  pass   :: print    => printQ
     procedure, public,  pass   :: defaults
     procedure, public,  pass   :: alloc
  end type qpProb

contains

  subroutine sqTitle( title )
    character(30), intent(out) :: title
  end subroutine sqTitle

  subroutine loadQ ( prob, ProbName, m, n, nnH, iObj, ObjAdd, &
                     neA, indA, locA, valA, bl, bu, ncObj, cObj, &
                     nNames, Names, hEtype, hs, x, pi, rc, usrHx )
    class(qpProb), target     :: prob
    character(8),  intent(in) :: ProbName
    integer(ip),   intent(in) :: m, n, nNames, ncObj, nnH, iObj, neA
    real(rp),      intent(in) :: ObjAdd
    character(8),  target     :: Names(:)
    integer(ip),   target     :: indA(:), locA(:), hEtype(:), hs(:)
    real(rp),      target     :: valA(:), bl(:), bu(:), cObj(:),               &
                                 x(:), pi(:), rc(:)
    procedure(iusrHx)         :: usrHx
  end subroutine loadQ

  subroutine loadQH ( prob, probName, m, n, nnH, iObj, ObjAdd, neA, indA,      &
                      locA, valA, bl, bu, ncObj, cObj, nNames, Names, hEtype,  &
                      hs, x, pi, rc, neH, indH, locH, valH )
    class(qpProb)            :: prob
    character(8), intent(in) :: probName
    integer(ip),  intent(in) :: m, n, nnH, ncObj, neA, iObj, nNames, neH
    real(rp),     intent(in) :: ObjAdd
    character(8), target     :: Names(:)
    integer(ip),  target     :: locA(:), indA(:), locH(:), indH(:),            &
                                hEtype(:), hs(:)
    real(rp),     target     :: cObj(:), bl(:), bu(:), valA(:), valH(:),       &
                                x(:), pi(:), rc(:)
  end subroutine loadQH

  subroutine solve ( prob, Start, INFO, nS, nInf, sInf, Obj )
    class(qpProb)                :: prob
    character*(*), intent(in)    :: Start
    integer(ip),   intent(inout) :: nInf
    integer(ip),   intent(out)   :: INFO, nS
    real(rp),      intent(out)   :: sInf, Obj
    INFO = 15
    Obj = 0.0D+0
  end subroutine solve

  subroutine defaults ( prob )
    class(qpProb) :: prob
    integer(ip) :: m, n, ne, ngObj, ngQP, nnH, iObj
  end subroutine defaults

  subroutine alloc ( prob, inform )
    class(qpProb)            :: prob
    integer(ip), intent(out) :: inform
  end subroutine alloc

  subroutine printQ ( prob )
    class(qpProb), target :: prob
  end subroutine printQ

end module sqModuleProb






