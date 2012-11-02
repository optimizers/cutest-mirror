! ( Last modified on 14 Jan 2001 at 19:03:33 )
      SUBROUTINE CISH ( N, X, IPROB, NNZH, LH, H, IRNH, ICNH )
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: N, IPROB, NNZH, LH
      INTEGER :: IRNH ( LH ), ICNH ( LH )
      REAL ( KIND = wp ) :: X ( N ), H( LH )

!  Compute the Hessian matrix of a specified problem function 
! (IPROB = 0 is the objective function, while IPROB > 0 is the 
!   IPROB-th constraint) of a problem initially written in 
!  Standard Input Format (SIF).

!  The upper triangle of the Hessian is stored in coordinate form,
!  i.e., the entry H(i) has row index IRNH(i) and column index ICNH(i)
!  for i = 1, ...., NNZH.

!  Based on the minimization subroutine LANCELOT/SBMIN
!  by Conn, Gould and Toint.

!  Nick Gould, for CGT productions,
!  August 1998.

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

!  integer variables from the PRFCTS common block.

      INTEGER :: NC2OF, NC2OG, NC2OH, NC2CF, NC2CG, NC2CH
      INTEGER :: NHVPR, PNC
      REAL :: SUTIME, STTIME
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
      COMMON / PRFCTS /  NC2OF, NC2OG, NC2OH, NC2CF, NC2CG, NC2CH, &
                         NHVPR, PNC, SUTIME, STTIME
      INTEGER :: IOUT
      COMMON / OUTPUT /  IOUT
      INTEGER :: NUMVAR, NUMCON
      COMMON / DIMS /    NUMVAR, NUMCON
      SAVE             / GLOBAL /, / LOCAL /, / CHARA /, / OUTPUT /, &
                       / DIMS   /, / PRFCTS /

!  Local variables

      INTEGER :: I, J, NCALCF, NCALCG, IFSTAT, IGSTAT
      INTEGER :: IG, LIH, LIWKH, LNXTRW, LINXTR, INFORM
      REAL ( KIND = wp ) :: ZERO, ONE, FTT
      PARAMETER ( ZERO = 0.0_wp, ONE = 1.0_wp )
      EXTERNAL :: RANGE 

!  Check input parameters.

      IF ( IPROB < 0 ) THEN
         WRITE( IOUT, 2020 )
         STOP
      END IF

!  Find group index IG of constraint IPROB.

      IF ( IPROB > 0 ) THEN
         IG = 0
         DO 10 I = 1, NG
            IF ( IWK( KNDOFC + I ) == IPROB ) THEN
              IG = I
              GO TO 20
            END IF
   10    CONTINUE
   20    CONTINUE
         IF ( IG == 0 ) THEN
            WRITE( IOUT, 2020 )
            STOP
         END IF
      END IF

!  Find which elements are involved in the required problem function.
!  Initialize the list to zero.

      DO 30 I = 1, NELNUM
        IWK( ICALCF + I ) = 0
   30 CONTINUE

!  If group IG is involved, record its elements.

      DO 50 IG = 1, NG
         IF ( IWK( KNDOFC + IG ) == IPROB ) THEN
            DO 40 J = IWK( ISTADG + IG ), IWK( ISTADG + IG + 1 ) - 1
               IWK( ICALCF + IWK( IELING + J ) ) = 1
   40       CONTINUE
         END IF
   50 CONTINUE

!  Only compute the elements which are involved in the required function.

      NCALCF = 0
      DO 80 I = 1, NELNUM
         IF ( IWK( ICALCF + I ) == 1 ) THEN
            NCALCF = NCALCF + 1
            IWK( ICALCF + NCALCF ) = I
         ELSE

!  If this is the first ever evaluation, initialize FUVALS.

            IF ( FIRSTG ) THEN
               FUVALS( I ) = ZERO
               DO 60 J = IWK( INTVAR + I ), IWK( INTVAR + I + 1 ) - 1
                  FUVALS( J ) = ZERO
   60          CONTINUE
               DO 70 J = IWK( ISTADH + I ), IWK( ISTADH + I + 1 ) - 1
                  FUVALS( J ) = ZERO
   70          CONTINUE
            END IF
         END IF
   80 CONTINUE

!  evaluate the element function values.

      CALL ELFUN ( FUVALS, X, WK( EPVALU + 1 ), NCALCF, &
                   IWK( ITYPEE + 1 ), IWK( ISTAEV + 1 ), &
                   IWK( IELVAR + 1 ), IWK( INTVAR + 1 ), &
                   IWK( ISTADH + 1 ), IWK( ISTEP + 1 ), &
                   IWK( ICALCF + 1 ),  &
                   LINTRE, LSTAEV, LELVAR, LNTVAR, LSTADH,  &
                   LNTVAR, LINTRE, LFUVAL, LVSCAL, LEPVLU,  &
                   1, IFSTAT )

!  evaluate the element function derivatives.

      CALL ELFUN ( FUVALS, X, WK( EPVALU + 1 ), NCALCF, &
                   IWK( ITYPEE + 1 ), IWK( ISTAEV + 1 ), &
                   IWK( IELVAR + 1 ), IWK( INTVAR + 1 ), &
                   IWK( ISTADH + 1 ), IWK( ISTEP + 1 ), &
                   IWK( ICALCF + 1 ),  &
                   LINTRE, LSTAEV, LELVAR, LNTVAR, LSTADH,  &
                   LNTVAR, LINTRE, LFUVAL, LVSCAL, LEPVLU,  &
                   3, IFSTAT )

!  Compute the list of groups involved in the required problem function.

      NCALCG = 0
      DO 130 IG = 1, NG
         IF ( IWK( KNDOFC + IG ) == IPROB ) THEN
            NCALCG = NCALCG + 1
            IWK( ICALCF + NCALCG ) = IG

!  compute the group argument values ft.

            FTT = - WK( B + IG )

!  include the contribution from the linear element.

            DO 110 J = IWK( ISTADA + IG ), IWK( ISTADA + IG + 1 ) - 1
               FTT = FTT + WK( A + J ) * X( IWK( ICNA + J ) )
  110       CONTINUE

!  include the contributions from the nonlinear elements.

            DO 120 J = IWK( ISTADG + IG ), IWK( ISTADG + IG + 1 ) - 1
               FTT = FTT + WK( ESCALE + J ) *  &
                      FUVALS( IWK( IELING + J ) )
  120          CONTINUE
            WK( FT + IG ) = FTT

!  Record the derivatives of trivial groups.

            IF ( LOGI( GXEQX + IG ) ) THEN
               WK( GVALS + NG + IG ) = ONE
               WK( GVALS + 2 * NG + IG ) = ZERO
            END IF
         ELSE

!  If this is the first ever evaluation, initialize GVALS.

            IF ( FIRSTG ) THEN
               WK( GVALS + NG + IG ) = ONE
               WK( GVALS + 2 * NG + IG ) = ZERO
            END IF
         END IF
  130 CONTINUE

!  Evaluate the group derivative values.

      IF ( .NOT. ALTRIV ) CALL GROUP ( WK ( GVALS + 1 ), NG, &
            WK( FT + 1 ), WK ( GPVALU + 1 ), NCALCG, &
            IWK( ITYPEG + 1 ), IWK( ISTGP + 1 ), &
            IWK( ICALCF + 1 ), &
            LCALCG, NG1, LCALCG, LCALCG, LGPVLU, &
            .TRUE., IGSTAT )

!  Define the real work space needed for ELGRD.
!  Ensure that there is sufficient space.

      IF ( LWK2 < NG ) THEN
         IF ( IOUT > 0 ) WRITE( IOUT, 2000 )
         STOP
      END IF
      IF ( NUMCON > 0 ) THEN

!  Change the group weightings to include the contributions from
!  the Lagrange multipliers.

         DO 140 IG = 1, NG
            I = IWK( KNDOFC + IG )
            IF ( I == IPROB ) THEN
               WK( LWKSTR + IG ) = WK( GSCALE + IG )
            ELSE
               WK( LWKSTR + IG ) = ZERO
            END IF
  140    CONTINUE

!  Compute the gradient value.

      CALL DELGRD( N, NG, FIRSTG, IWK( ICNA + 1 ), LICNA, &
                      IWK( ISTADA + 1 ), LSTADA, IWK( IELING + 1 ), &
                      LELING, IWK( ISTADG + 1 ), LSTADG, &
                      IWK( ITYPEE + 1 ), LINTRE, &
                      IWK( ISTAEV + 1 ), LSTAEV, IWK( IELVAR + 1 ), &
                      LELVAR, IWK( INTVAR + 1 ), LNTVAR, &
                      IWK( LSVGRP + 1 ), &
                      LNVGRP, IWK( LSTAJC + 1 ), LNSTJC, &
                      IWK( LSTAGV + 1 ), LNSTGV, WK( A + 1 ), LA, &
                      WK( GVALS + NG + 1 ), LGVALS, &
                      FUVALS, LNGUVL, FUVALS( LGGFX + 1 ), &
                      WK( LWKSTR + 1 ), NG, &
                      WK( ESCALE + 1 ), LESCAL, FUVALS( LGRJAC + 1 ), &
                      LNGRJC, WK( WRK + 1 ), WK( WRK + N + 1 ), MAXSEL, &
                      LOGI( GXEQX + 1 ), LGXEQX, &
                      LOGI( INTREP + 1 ), LINTRE, RANGE )
      ELSE

!  Compute the gradient value.

      CALL DELGRD( N, NG, FIRSTG, IWK( ICNA + 1 ), LICNA, &
                      IWK( ISTADA + 1 ), LSTADA, IWK( IELING + 1 ), &
                      LELING, IWK( ISTADG + 1 ), LSTADG, &
                      IWK( ITYPEE + 1 ), LINTRE, &
                      IWK( ISTAEV + 1 ), LSTAEV, IWK( IELVAR + 1 ), &
                      LELVAR, IWK( INTVAR + 1 ), LNTVAR, &
                      IWK( LSVGRP + 1 ), &
                      LNVGRP, IWK( LSTAJC + 1 ), LNSTJC, &
                      IWK( LSTAGV + 1 ), LNSTGV, WK( A + 1 ), LA, &
                      WK( GVALS + NG + 1 ), LGVALS, &
                      FUVALS, LNGUVL, FUVALS( LGGFX + 1 ), &
                      WK( GSCALE + 1 ), LGSCAL, &
                      WK( ESCALE + 1 ), LESCAL, FUVALS( LGRJAC + 1 ), &
                      LNGRJC, WK( WRK + 1 ), WK( WRK + N + 1 ), MAXSEL, &
                      LOGI( GXEQX + 1 ), LGXEQX, &
                      LOGI( INTREP + 1 ), LINTRE, RANGE )
      END IF
      FIRSTG = .FALSE.

!  Define the real work space needed for ASMBL.
!  Ensure that there is sufficient space.

      IF ( NUMCON > 0 ) THEN
         IF ( LWK2 < N + 3 * MAXSEL + NG ) THEN
            IF ( IOUT > 0 ) WRITE( IOUT, 2000 )
            STOP
         END IF
      ELSE
         IF ( LWK2 < N + 3 * MAXSEL ) THEN
            IF ( IOUT > 0 ) WRITE( IOUT, 2000 )
            STOP
         END IF
      END IF

!  Define the integer work space needed for ASMBL.
!  Ensure that there is sufficient space.

      LIWKH = LIWK2 - N
      LINXTR = LIWKH / 2
      IF ( LINXTR < N ) THEN
         IF ( IOUT > 0 ) WRITE( IOUT, 2010 )
         STOP
      END IF

!  Set starting addresses for partitions of the integer workspace.

      LIH = LH
      LNXTRW = 0
      DO 190 I = 1, N
         IWK( IVAR + I ) = I
  190 CONTINUE

!  Assemble the Hessian.

      IF ( NUMCON > 0 ) THEN
      CALL DASMBL( N, NG, MAXSEL, N, LH, LIH, NNZH, &
                   N, IWK( IVAR + 1), IWK( ISTADH + 1 ), LSTADH, &
                   IWK( ICNA + 1 ), LICNA, &
                   IWK( ISTADA + 1 ), LSTADA, IWK( INTVAR + 1 ), LNTVAR, &
                   IWK( IELVAR + 1 ), LELVAR, IWK( IELING + 1 ), LELING, &
                   IWK( ISTADG + 1 ), LSTADG, IWK( ISTAEV + 1 ), LSTAEV, &
                   IWK( LSTAGV + 1 ), LNSTGV, IWK( LSVGRP + 1 ), LNVGRP, &
                   IRNH, ICNH, IWK( LSEND + LNXTRW + 1 ), LINXTR, &
                   IWK( LSEND + LIWKH + 1 ), N, &
                   WK( A + 1 ), LA, FUVALS, LNGUVL, FUVALS, LNHUVL, &
                   WK( GVALS + NG + 1 ), WK( GVALS + 2 * NG + 1 ), &
                   WK( LWKSTR + 1 ), WK( ESCALE + 1 ), LESCAL, &
                   H, WK ( LWKSTR + NG + 1 ), LWK2 - NG, &
                   LOGI( GXEQX + 1 ), LGXEQX, LOGI( INTREP + 1 ), &
                   LINTRE, IWK( ITYPEE + 1 ), LINTRE, &
                   RANGE, 1, IOUT, .FALSE., I, INFORM,  &
                   .FALSE., .FALSE. )
      ELSE
      CALL DASMBL( N, NG, MAXSEL, N, LH, LIH, NNZH, &
                   N, IWK( IVAR + 1), IWK( ISTADH + 1 ), LSTADH, &
                   IWK( ICNA + 1 ), LICNA, &
                   IWK( ISTADA + 1 ), LSTADA, IWK( INTVAR + 1 ), LNTVAR, &
                   IWK( IELVAR + 1 ), LELVAR, IWK( IELING + 1 ), LELING, &
                   IWK( ISTADG + 1 ), LSTADG, IWK( ISTAEV + 1 ), LSTAEV, &
                   IWK( LSTAGV + 1 ), LNSTGV, IWK( LSVGRP + 1 ), LNVGRP, &
                   IRNH, ICNH, IWK( LSEND + LNXTRW + 1 ), LINXTR, &
                   IWK( LSEND + LIWKH + 1 ), N, &
                   WK( A + 1 ), LA, FUVALS, LNGUVL, FUVALS, LNHUVL, &
                   WK( GVALS + NG + 1 ), WK( GVALS + 2 * NG + 1 ), &
                   WK( GSCALE + 1 ), WK( ESCALE + 1 ), LESCAL, &
                   H, WK ( LWKSTR + 1 ), LWK2 - NG, &
                   LOGI( GXEQX + 1 ), LGXEQX, LOGI( INTREP + 1 ), &
                   LINTRE, IWK( ITYPEE + 1 ), LINTRE, &
                   RANGE, 1, IOUT, .FALSE., I, INFORM, &
                   .FALSE., .FALSE. )
      END IF

!  Check that there is sufficient integer workspace.

      IF ( INFORM > 0 ) THEN
         IF ( IOUT > 0 ) WRITE( IOUT, 2010 )
         STOP
      END IF

!  Update the counters for the report tool.

      IF( IPROB == 0 ) THEN
         NC2OH = NC2OH + 1
      ELSE 
         NC2CH = NC2CH + 1
      ENDIF 
      RETURN

! Non-executable statements.

 2000 FORMAT( ' ** SUBROUTINE CISH: Increase the size of WK ' )
 2010 FORMAT( ' ** SUBROUTINE CISH: Increase the sizes of IWK and LH' )
 2020 FORMAT( ' ** SUBROUTINE CISH: invalid problem index IPROB ' )

!  end of CISH.

      END
