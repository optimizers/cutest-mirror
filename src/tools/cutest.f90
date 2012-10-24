! THIS VERSION: CUTEST 1.0 - 24/10/2012 AT 13:00 GMT.

!-*-*-*-*-*-*-*-*-*-*-*- S I F D E C O D E   M O D U L E -*-*-*-*-*-*-*-*-*-*-

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

!  initial integer workspace

      INTEGER, PARAMETER :: liwk = 40000000

!  initial real precision workspace

      INTEGER, PARAMETER :: lwk = 50000000

!  initial logical workspace

      INTEGER, PARAMETER :: llogic = 2000000

!  initial character workspace

      INTEGER, PARAMETER :: lchara = 2000000

!  initial length of workspace to store the problem's function and 
!  derivatives values

      INTEGER,  PARAMETER :: lfuval = 8000000

!-------------------------------------------------
!  D e r i v e d   t y p e   d e f i n i t i o n s
!-------------------------------------------------

!  - - - - - - - - - -
!   data derived type
!  - - - - - - - - - -

      TYPE, PUBLIC :: CUTEST_data_type

!  Integer variables from the GLOBAL common block.

        INTEGER :: ng, nelnum, ngel, nvars, nnza, ngpvlu
        INTEGER :: nepvlu, ng1, nel1, istadg, istgp, istada
        INTEGER :: istaev, istep, itypeg, kndofc, itypee
        INTEGER :: ieling, ielvar, icna, istadh, intvar, ivar
        INTEGER :: icalcf, itypev, iwrk, a, b
        INTEGER :: u, gpvalu, epvalu
        INTEGER :: escale, gscale, vscale, gvals, xt, dgrad
        INTEGER :: q, wrk, intrep, gxeqx, gnames, vnames
        INTEGER :: lo, ch, liwork, lwork, ngng, ft
        INTEGER :: la, lb, nobjgr, lu, lelvar
        INTEGER :: lstaev, lstadh, lntvar, lcalcf
        INTEGER :: leling, lintre, lft, lgxeqx, lstadg, lgvals
        INTEGER :: licna, lstada, lkndof, lgpvlu, lepvlu
        INTEGER :: lgscal, lescal, lvscal, lcalcg            

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

!  workspace arrays

        INTEGER :: IWK( liwk )
        LOGICAL :: LOGI( llogic )
        CHARACTER ( len = 10 ):: CHA( lchara )
        REAL ( KIND = wp ):: WK lwk )
        REAL ( KIND = wp ) :: FUVALS( lfuval )

      END TYPE CUTEST_data_type

    END MODULE CUTEST
