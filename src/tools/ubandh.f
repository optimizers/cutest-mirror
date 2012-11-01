! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE UBANDH( N, GOTH, X, NSEMIB, BANDH, LBANDH )
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: N, NSEMIB, LBANDH
      LOGICAL :: GOTH
      REAL ( KIND = wp ) :: X( N ), BANDH( 0: LBANDH, N )

!  Compute the portion of the Hessian matrix of a group partially
!  separable function which lies within a band of given semi-bandwidth.
!  The diagonal and subdiagonal entries in column i are stored in
!  locations BAND( j, i ), where j = 0 for the diagonal and j > 0 for
!  the j-th subdiagonal entry, j = 1, ..., MIN( NSEMIB, N - i ) and
!  NSEMIB is the requested semi-bandwidth.

!  Based on the minimization subroutine LANCELOT/SBMIN
!  by Conn, Gould and Toint.

!  Nick Gould, for CGT productions.
!  August 1993.

      INTEGER :: LIWK, LWK, LFUVAL, LLOGIC, LCHARA

! ---------------------------------------------------------------------

!  Parameters whose value might be changed by the user:

!  The following parameters define the sizes of problem
!  dependent arrays. These may be changed by the user to
!  suit a particular problem or system configuration.

!  The TOOLS will issue error messages if any of these sizes
!  is too small, telling which parameter to increase.

! ---------------------------------------------------------------------

      INCLUDE 'tools.siz'

      INTEGER :: IWK( LIWK )
      LOGICAL :: LOGI ( LLOGIC )
      CHARACTER ( LEN = 10 ) :: CHA ( LCHARA )
      REAL ( KIND = wp ) :: WK ( LWK )
      REAL ( KIND = wp ) :: FUVALS ( LFUVAL )

! ---------------------------------------------------------------------

!  End of parameters which might be changed by the user.

! ---------------------------------------------------------------------

!  integer variables from the GLOBAL common block.

      INTEGER :: NG, NELNUM, NGEL, NVARS, NNZA, NGPVLU
      INTEGER :: NEPVLU, NG1, NEL1, ISTADG, ISTGP, ISTADA
      INTEGER :: ISTAEV, ISTEP, ITYPEG, KNDOFC, ITYPEE
      INTEGER :: IELING, IELVAR, ICNA, ISTADH, INTVAR, IVAR
      INTEGER :: ICALCF, ITYPEV, IWRK, A, B
      INTEGER :: U, GPVALU, EPVALU
      INTEGER :: ESCALE, GSCALE, VSCALE, GVALS, XT, DGRAD
      INTEGER :: Q, WRK, INTREP, GXEQX, GNAMES, VNAMES
      INTEGER :: LO, CH, LIWORK, LWORK, NGNG, FT
      INTEGER :: LA, LB, NOBJGR, LU, LELVAR
      INTEGER :: LSTAEV, LSTADH, LNTVAR, LCALCF
      INTEGER :: LELING, LINTRE, LFT, LGXEQX, LSTADG, LGVALS
      INTEGER :: LICNA, LSTADA, LKNDOF, LGPVLU, LEPVLU
      INTEGER :: LGSCAL, LESCAL, LVSCAL, LCALCG

!  integer variables from the LOCAL common block.

      INTEGER :: LFXI, LGXI, LHXI, LGGFX, LDX, LGRJAC
      INTEGER :: LQGRAD, LBREAK, LP, LXCP, LX0, LGX0
      INTEGER :: LDELTX, LBND, LWKSTR, LSPTRS, LSELTS, LINDEX
      INTEGER :: LSWKSP, LSTAGV, LSTAJC, LIUSED, LFREEC
      INTEGER :: LNNONZ, LNONZ2, LSYMMD, LSYMMH
      INTEGER :: LSLGRP, LSVGRP, LGCOLJ, LVALJR, LSEND
      INTEGER :: LNPTRS, LNELTS, LNNDEX, LNWKSP, LNSTGV
      INTEGER :: LNSTJC, LNIUSE, LNFREC, LNNNON, LNNNO2, LNSYMD
      INTEGER :: LNSYMH, LNLGRP, LNVGRP, LNGCLJ, LNVLJR, LNQGRD
      INTEGER :: LNBRAK, LNP, LNBND, LNFXI, LNGXI, LNGUVL
      INTEGER :: LNHXI, LNHUVL, LNGGFX, LNDX, LNGRJC, LIWK2
      INTEGER :: LWK2, MAXSIN, NINVAR, MAXSEL
      INTEGER :: NTYPE, NSETS, LSTYPE, LSSWTR, LSSIWT, LSIWTR
      INTEGER :: LSWTRA, LNTYPE, LNSWTR, LNSIWT, LNIWTR
      INTEGER :: LNWTRA, LSISET, LSSVSE, LNISET, LNSVSE
      LOGICAL :: ALTRIV, FIRSTG
      COMMON / GLOBAL /  IWK, WK, FUVALS, LOGI, &
                         NG, NELNUM, NGEL, NVARS, NNZA, NGPVLU, &
                         NEPVLU, NG1, NEL1, ISTADG, ISTGP, ISTADA, &
                         ISTAEV, ISTEP, ITYPEG, KNDOFC, ITYPEE, &
                         IELING, IELVAR, ICNA, ISTADH, INTVAR, IVAR, &
                         ICALCF, ITYPEV, IWRK, A, B, &
                         U, GPVALU, EPVALU, &
                         ESCALE, GSCALE, VSCALE, GVALS, XT, DGRAD, &
                         Q, WRK, INTREP, GXEQX, GNAMES, VNAMES, &
                         LO, CH, LIWORK, LWORK, NGNG, FT, &
                         ALTRIV, FIRSTG, &
                         LA, LB, NOBJGR, LU, LELVAR, &
                         LSTAEV, LSTADH, LNTVAR, LCALCF, &
                         LELING, LINTRE, LFT, LGXEQX, LSTADG, LGVALS, &
                         LICNA, LSTADA, LKNDOF, LGPVLU, LEPVLU, &
                         LGSCAL, LESCAL, LVSCAL, LCALCG
      COMMON / CHARA /   CHA
      COMMON / LOCAL /   LFXI, LGXI, LHXI, LGGFX, LDX, LGRJAC, &
                         LQGRAD, LBREAK, LP, LXCP, LX0, LGX0, &
                         LDELTX, LBND, LWKSTR, LSPTRS, LSELTS, LINDEX, &
                         LSWKSP, LSTAGV, LSTAJC, LIUSED, LFREEC, &
                         LNNONZ, LNONZ2, LSYMMD, LSYMMH, &
                         LSLGRP, LSVGRP, LGCOLJ, LVALJR, LSEND, &
                         LNPTRS, LNELTS, LNNDEX, LNWKSP, LNSTGV, &
                         LNSTJC, LNIUSE, LNFREC, LNNNON, LNNNO2, LNSYMD, &
                         LNSYMH, LNLGRP, LNVGRP, LNGCLJ, LNVLJR, LNQGRD, &
                         LNBRAK, LNP, LNBND, LNFXI, LNGXI, LNGUVL, &
                         LNHXI, LNHUVL, LNGGFX, LNDX, LNGRJC, LIWK2, &
                         LWK2, MAXSIN, NINVAR, MAXSEL, NTYPE, &
                         NSETS, LSTYPE, LSSWTR, LSSIWT, LSIWTR, &
                         LSWTRA, LNTYPE, LNSWTR, LNSIWT, LNIWTR, &
                         LNWTRA, LSISET, LSSVSE, LNISET, LNSVSE
      INTEGER :: IOUT
      COMMON / OUTPUT /  IOUT
      SAVE             / GLOBAL /, / LOCAL /, / CHARA /, / OUTPUT /

!  Local variables

      INTEGER :: I, IG, J, LH2, NSEMIW, INFORM, NNZH
      INTEGER :: IFSTAT, IGSTAT
      INTEGER :: IRNH(   1 ), ICNH(   1 )
      REAL ( KIND = wp ) :: ZERO, ONE, FTT
      PARAMETER ( ZERO = 0.0_wp, ONE = 1.0_wp )
      EXTERNAL :: RANGE 

!  Check that there is room for the band matrix.

      NSEMIW = MAX( 0, MIN( NSEMIB, N - 1 ) )
      IF ( LBANDH < NSEMIW ) THEN
         IF ( IOUT > 0 ) WRITE( IOUT, 2020 )
         STOP
      END IF

!  there are non-trivial group functions.

      DO 10 I = 1, MAX( NELNUM, NG )
        IWK( ICALCF + I ) = I
   10 CONTINUE

!  evaluate the element function values.

      CALL ELFUN ( FUVALS, X, WK( EPVALU + 1 ), NELNUM, &
                   IWK( ITYPEE + 1 ), IWK( ISTAEV + 1 ), &
                   IWK( IELVAR + 1 ), IWK( INTVAR + 1 ), &
                   IWK( ISTADH + 1 ), IWK( ISTEP + 1 ), &
                   IWK( ICALCF + 1 ),  &
                   LINTRE, LSTAEV, LELVAR, LNTVAR, LSTADH,  &
                   LNTVAR, LINTRE, LFUVAL, LVSCAL, LEPVLU,  &
                   1, IFSTAT )

!  evaluate the element function values.

      CALL ELFUN ( FUVALS, X, WK( EPVALU + 1 ), NELNUM, &
                   IWK( ITYPEE + 1 ), IWK( ISTAEV + 1 ), &
                   IWK( IELVAR + 1 ), IWK( INTVAR + 1 ), &
                   IWK( ISTADH + 1 ), IWK( ISTEP + 1 ), &
                   IWK( ICALCF + 1 ),  &
                   LINTRE, LSTAEV, LELVAR, LNTVAR, LSTADH,  &
                   LNTVAR, LINTRE, LFUVAL, LVSCAL, LEPVLU,  &
                   3, IFSTAT )

!  compute the group argument values ft.

      DO 90 IG = 1, NG
         FTT = - WK( B + IG )

!  include the contribution from the linear element.

         DO 30 J = IWK( ISTADA + IG ), IWK( ISTADA + IG + 1 ) - 1
            FTT = FTT + WK( A + J ) * X( IWK( ICNA + J ) )
   30    CONTINUE

!  include the contributions from the nonlinear elements.

         DO 60 J = IWK( ISTADG + IG ), IWK( ISTADG + IG + 1 ) - 1
            FTT = FTT + WK( ESCALE + J ) * FUVALS( IWK( IELING + J ) )
   60    CONTINUE
         WK( FT + IG ) = FTT

!  Record the derivatives of trivial groups.

         IF ( LOGI( GXEQX + IG ) ) THEN
            WK( GVALS + NG + IG ) = ONE
            WK( GVALS + 2 * NG + IG ) = ZERO
         END IF
   90 CONTINUE

!  evaluate the group derivative values.

      IF ( .NOT. ALTRIV ) CALL GROUP ( WK ( GVALS + 1 ), NG, &
            WK( FT + 1 ), WK ( GPVALU + 1 ), NG, &
            IWK( ITYPEG + 1 ), IWK( ISTGP + 1 ), &
            IWK( ICALCF + 1 ), &
            LCALCG, NG1, LCALCG, LCALCG, LGPVLU, &
            .TRUE., IGSTAT )

!  Compute the gradient value.

      CALL DELGRD( N, NG, FIRSTG, IWK( ICNA + 1 ), LICNA, &
                   IWK( ISTADA + 1 ), LSTADA, IWK( IELING + 1 ), &
                   LELING, IWK( ISTADG + 1 ), LSTADG, &
                   IWK( ITYPEE + 1 ), LINTRE, &
                   IWK( ISTAEV + 1 ), LSTAEV, IWK( IELVAR + 1 ), &
                   LELVAR, IWK( INTVAR + 1 ), LNTVAR, IWK( LSVGRP + 1 ), &
                   LNVGRP, IWK( LSTAJC + 1 ), LNSTJC, &
                   IWK( LSTAGV + 1 ), LNSTGV, WK( A + 1 ), LA, &
                   WK( GVALS + NG + 1 ), LGVALS, &
                   FUVALS, LNGUVL, FUVALS( LGGFX + 1 ), &
                   WK( GSCALE + 1 ), LGSCAL, &
                   WK( ESCALE + 1 ), LESCAL, FUVALS( LGRJAC + 1 ), &
                   LNGRJC, WK( WRK + 1 ), WK( WRK + N + 1 ), MAXSEL, &
                   LOGI( GXEQX + 1 ), LGXEQX, &
                   LOGI( INTREP + 1 ), LINTRE, RANGE )
      FIRSTG = .FALSE.

!  Define the real work space needed for ASMBL.
!  Ensure that there is sufficient space.

      LH2 = N * ( NSEMIW + 1 )
      IF ( LWK2 < LH2 + N + 3 * MAXSEL ) THEN
         IF ( IOUT > 0 ) WRITE( IOUT, 2000 )
         STOP
      END IF

!  Set starting addresses for partitions of the integer workspace.

      DO 100 I = 1, N
         IWK( IVAR + I ) = I
  100 CONTINUE

!  Assemble the Hessian.

      CALL DASMBL( N, NG, MAXSEL, NSEMIW, LH2, 1, NNZH, &
                   N, IWK( IVAR + 1), IWK( ISTADH + 1 ), LSTADH, &
                   IWK( ICNA + 1 ), LICNA, &
                   IWK( ISTADA + 1 ), LSTADA, IWK( INTVAR + 1 ), LNTVAR, &
                   IWK( IELVAR + 1 ), LELVAR, IWK( IELING + 1 ), LELING, &
                   IWK( ISTADG + 1 ), LSTADG, IWK( ISTAEV + 1 ), LSTAEV, &
                   IWK( LSTAGV + 1 ), LNSTGV, IWK( LSVGRP + 1 ), LNVGRP, &
                   IRNH, ICNH, IWK( LSEND + 1 ), 1, &
                   IWK( LIWK2 - N + 1 ), N, &
                   WK( A + 1 ), LA, FUVALS, LNGUVL, FUVALS, LNHUVL, &
                   WK( GVALS + NG + 1 ), WK( GVALS + 2 * NG + 1 ), &
                   WK( GSCALE + 1 ), WK( ESCALE + 1 ), LESCAL, &
                   WK( LWKSTR + 1 ), WK ( LWKSTR + LH2 + 1 ), &
                   LWK2 - LH2, LOGI( GXEQX + 1 ), &
                   LGXEQX, LOGI( INTREP + 1 ), LINTRE, &
                   IWK( ITYPEE + 1 ), LINTRE, &
                   RANGE, 1, IOUT, .TRUE., I, INFORM, .FALSE., &
                   .FALSE. )

!  Check that there is sufficient integer workspace.

      IF ( INFORM > 0 ) THEN
         IF ( IOUT > 0 ) WRITE( IOUT, 2010 )
         STOP
      END IF

!  Place the diagonal entries in their correct positions.

      DO 110 J = 1, N
         BANDH( 0, J ) = WK( LWKSTR + J )
  110 CONTINUE

!  Now, place the subdiagonal entries in their correct positions.

      LH2 = LWKSTR + N
      DO 130 J = 1, N
         DO 120 I = 1, MIN( NSEMIW, N - J )
            BANDH( I, J ) = WK( LH2 + I )
  120    CONTINUE
         LH2 = LH2 + NSEMIW
  130 CONTINUE
      RETURN

! Non-executable statements.

 2000 FORMAT( ' ** SUBROUTINE UBANDH: Increase the size of WK ' )
 2010 FORMAT( ' ** SUBROUTINE UBANDH: Increase the size of IWK ' )
 2020 FORMAT( ' ** SUBROUTINE UBANDH: array dimension LBANDH should be', &
              ' larger than NSEMIB' )

!  end of UBANDH.

      END

