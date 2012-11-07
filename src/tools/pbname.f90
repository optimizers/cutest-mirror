! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE PBNAME( n, PNAME )
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n
      CHARACTER ( LEN = 10 ) :: PNAME

!  Obtain the name of the problem

!  D. Orban, for GOT productions, Aug 2005, from the version of
!  Nick Gould, for CGT productions.
!  September 1992.

      INTEGER :: liwk, lwk, lfuval, llogic, lchara

! ---------------------------------------------------------------------

!  Parameters whose value might be changed by the user:

!  The following parameters define the sizes of problem
!  dependent arrays. These may be changed by the user to
!  suit a particular problem or system configuration.

!  The TOOLS will issue error messages if any of these sizes
!  is too small, telling which parameter to increase.

! ---------------------------------------------------------------------

      INCLUDE 'tools.siz'

      INTEGER :: IWK( liwk )
      LOGICAL :: LOGI ( llogic )
      CHARACTER ( LEN = 10 ) :: CHA ( lchara )
      REAL ( KIND = wp ) :: WK ( lwk )
      REAL ( KIND = wp ) :: FUVALS ( lfuval )

! ---------------------------------------------------------------------

!  End of parameters which might be changed by the user.

! ---------------------------------------------------------------------

!  integer variables from the GLOBAL common block.

      INTEGER :: ng, nel, ntotel, nvrels, nnza, ngpvlu
      INTEGER :: nepvlu, ng1, nel1, ISTADG, ISTGP, ISTADA
      INTEGER :: ISTAEV, ISTEP, ITYPEG, KNDOFC, ITYPEE
      INTEGER :: IELING, IELVAR, ICNA, ISTADH, INTVAR, IVAR
      INTEGER :: ICALCF, ITYPEV, IWRK, A, B
      INTEGER :: U, GPVALU, EPVALU
      INTEGER :: ESCALE, GSCALE, VSCALE, GVALS, XT, DGRAD
      INTEGER :: Q, WRK, intrep, gxeqx, GNAMES, VNAMES
      INTEGER :: lo, ch, liwork, lwork, ngng, FT
      INTEGER :: la, lb, nobjgr, lu, lelvar
      INTEGER :: lstaev, lstadh, lntvar, lcalcf
      INTEGER :: leling, lintre, lft, lgxeqx, lstadg, lgvals
      INTEGER :: licna, lstada, lkndof, lgpvlu, lepvlu
      INTEGER :: lgscal, lescal, lvscal, lcalcg

!  integer variables from the LOCAL common block.

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
      INTEGER :: lwk2, maxsin, ninvar, maxsel
      INTEGER :: ntype, nsets, lstype, lsswtr, lssiwt, lsiwtr
      INTEGER :: lswtra, lntype, lnswtr, lnsiwt, lniwtr
      INTEGER :: lnwtra, lsiset, lssvse, lniset, lnsvse
      LOGICAL :: ALTRIV, FIRSTG
      COMMON / GLOBAL /  IWK, WK, FUVALS, LOGI, &
                         ng, nel, ntotel, nvrels, nnza, ngpvlu, &
                         nepvlu, ng1, nel1, ISTADG, ISTGP, ISTADA, &
                         ISTAEV, ISTEP, ITYPEG, KNDOFC, ITYPEE, &
                         IELING, IELVAR, ICNA, ISTADH, INTVAR, IVAR, &
                         ICALCF, ITYPEV, IWRK, A, B, &
                         U, GPVALU, EPVALU, &
                         ESCALE, GSCALE, VSCALE, GVALS, XT, DGRAD, &
                         Q, WRK, intrep, gxeqx, GNAMES, VNAMES, &
                         lo, ch, liwork, lwork, ngng, FT, &
                         ALTRIV, FIRSTG, &
                         la, lb, nobjgr, lu, lelvar, &
                         lstaev, lstadh, lntvar, lcalcf, &
                         leling, lintre, lft, lgxeqx, lstadg, lgvals, &
                         licna, lstada, lkndof, lgpvlu, lepvlu, &
                         lgscal, lescal, lvscal, lcalcg
      COMMON / CHARA /   CHA
      COMMON / LOCAL /   lfxi, lgxi, lhxi, lggfx, ldx, lgrjac, &
                         lqgrad, lbreak, lp, lxcp, lx0, lgx0, &
                         ldeltx, lbnd, lwkstr, lsptrs, lselts, lindex, &
                         lswksp, lstagv, lstajc, liused, lfreec, &
                         lnnonz, lnonz2, lsymmd, lsymmh, &
                         lslgrp, lsvgrp, lgcolj, lvaljr, lsend, &
                         lnptrs, lnelts, lnndex, lnwksp, lnstgv, &
                         lnstjc, lniuse, lnfrec, lnnnon, lnnno2, lnsymd, &
                         lnsymh, lnlgrp, lnvgrp, lngclj, lnvljr, lnqgrd, &
                         lnbrak, lnp, lnbnd, lnfxi, lngxi, lnguvl, &
                         lnhxi, lnhuvl, lnggfx, lndx, lngrjc, liwk2, &
                         lwk2, maxsin, ninvar, maxsel, ntype, &
                         nsets, lstype, lsswtr, lssiwt, lsiwtr, &
                         lswtra, lntype, lnswtr, lnsiwt, lniwtr, &
                         lnwtra, lsiset, lssvse, lniset, lnsvse
      SAVE             / GLOBAL /, / LOCAL /, / CHARA /

!  local variables.

      INTEGER :: i

!  Set the problem name.

      PNAME = CHA( VNAMES + n + 1 )
      RETURN

!  end of PBNAME.

      END
