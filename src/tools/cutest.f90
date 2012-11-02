! THIS VERSION: CUTEST 1.0 - 24/10/2012 AT 13:00 GMT.

!-*-*-*-*-*-*-*-*-*-*-*-*-*- C U T E S T   M O D U L E -*-*-*-*-*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 77 version originally released as part of CUTE, December 1990
!   Became separate package CUTEr, April 2004
!   Updated fortran 2003 version released October 2012

!  For full documentation, see 
!   http://galahad.rl.ac.uk/galahad-www/specs.html

    MODULE CUTEST

      IMPLICIT NONE

      PRIVATE

!--------------------
!   P r e c i s i o n
!--------------------

      INTEGER, PARAMETER :: sp = KIND( 1.0E+0 )
      INTEGER, PARAMETER :: dp = KIND( 1.0D+0 )
      INTEGER, PARAMETER :: wp = dp

!----------------------
!   P a r a m e t e r s
!----------------------

      REAL ( KIND = wp ), PARAMETER :: zero = 0.0_wp
      REAL ( KIND = wp ), PARAMETER :: one = 1.0_wp

!-------------------------------------------------
!  D e r i v e d   t y p e   d e f i n i t i o n s
!-------------------------------------------------

!  - - - - - - - - - -
!   data derived type
!  - - - - - - - - - -

      TYPE, PUBLIC :: CUTEST_data_type

!  Integer variables from the GLOBAL common block.

        INTEGER :: ng, nelnum, ngel, nvars, nnza, ngpvlu, nepvlu, ng1, nel1
        INTEGER :: lo, ch, liwork, lwork, ngng, la, lb, nobjgr, lu, lelvar
        INTEGER :: lstaev, lstadh, lntvar, lcalcf, leling, lintre, lft
        INTEGER :: lgxeqx, licna, lstada, lkndof, lgpvlu, lepvlu
        INTEGER :: lstadg, lgvals, lgscal, lescal, lvscal, lcalcg            
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: ISTADG
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: ISTGP
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: ISTADA
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: ISTAEV
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: ISTEP
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: ITYPEG
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: KNDOFC
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: ITYPEE
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: IELING
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: IELVAR
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: ICNA
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: ISTADH
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: INTVAR
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: IVAR
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: ICALCF
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: ITYPEV
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: IWORK
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: A
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: B
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: U
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: GPVALU
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: EPVALU
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: ESCALE
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: GSCALE
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: VSCALE
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: GVALS
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: XT
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DGRAD
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: Q
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: FT
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: WRK
        LOGICAL, ALLOCATABLE, DIMENSION( : ) :: INTREP
        LOGICAL, ALLOCATABLE, DIMENSION( : ) :: GXEQX
        LOGICAL, ALLOCATABLE, DIMENSION( : ) :: LOGIC
        CHARACTER ( LEN = 10 ), ALLOCATABLE, DIMENSION( : ) :: GNAMES
        CHARACTER ( LEN = 10 ), ALLOCATABLE, DIMENSION( : ) :: VNAMES
        CHARACTER ( LEN = 10 ) :: pname
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: FUVALS
!       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: FXI
!       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: GXI
!       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: HXI
!       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: GGFX
!       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DX
!       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: GRJAC

!  Integer variables from the LOCAL common block.

        INTEGER :: lfxi, lgxi, lhxi, lggfx, ldx, lgrjac
        INTEGER :: lqgrad, lbreak, lp, lxcp, lx0, lgx0
        INTEGER :: ldeltx, lbnd, lwkstr, lsptrs, lselts, lindex
        INTEGER :: lswksp, lstagv, lstajc, liused, lfreec
        INTEGER :: lnnonz, lnonz2, lsymmd, lsymmh
        INTEGER :: lslgrp, lsvgrp, lgcolj, lvaljr, lsend
        INTEGER :: lnptrs, lnelts, lnndex, lnwksp, lnstgv
        INTEGER :: lnstjc, lniuse, lnfrec, lnnnon, lnnno2, lnsymd
        INTEGER :: lnsymh, lnlgrp, lnvgrp, lngclj, lnvljr, lnqgrd
        INTEGER :: lnbrak, lnp, lnbnd, lnfxi, lngxi, lnguvl
        INTEGER :: lnhxi, lnhuvl, lnggfx, lndx, lngrjc, liwk2
        INTEGER :: lwk2, maxsin, ninvar, maxsel, lniwtr
        INTEGER :: ntype, nsets, lstype, lsswtr, lssiwt
        INTEGER :: lsiwtr, lswtra, lntype, lnswtr, lnsiwt
        INTEGER :: lnwtra, lsiset, lssvse, lniset, lnsvse
        LOGICAL :: altriv, firstg                                 

!  variables from the PRFCTS common block

        INTEGER :: nc2of, nc2og, nc2oh, nc2cf, nc2cg, nc2ch, nhvpr, pnc
        REAL :: sutime, sttime

!  variables from the NNVARS common block

        INTEGER :: nnov, nnjv 

!  variables from the OUTPUT common block

        INTEGER :: iout2

!  variables from the DIMS common block

        INTEGER :: numvar, numcon

      END TYPE CUTEST_data_type

    END MODULE CUTEST
