! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE CDIMSH( NNZH )
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: NNZH

!  Compute the space required to store the Hessian matrix of the 
!  Lagrangian function of a problem initially written in 
!  Standard Input Format (SIF).

!  The upper triangle of the Hessian is stored in coordinate form,
!  i.e., the entry H(i) has row index IRNH(i)
!  for i = 1, ...., NNZH.

!  Based on the minimization subroutine LANCELOT/SBMIN
!  by Conn, Gould and Toint.

!  Nick Gould, for CGT productions,
!  August 1999.

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

      INTEGER :: NXTRW1, NXTRW2, LNXTRW, LIRNH, IRNH
      INTEGER :: I, II,  IG, J,  JJ, K,  L, IEL, IELL, INEXT 
      INTEGER :: NEWPT,  IG1,    LISTVS, LISTVE, ISTART

!  Define the integer work space needed for ASMBI.
!  Ensure that there is sufficient space

      NNZH = 0
      NEWPT = NUMVAR + 1
      LIRNH = ( LIWK2 - 2 * NUMVAR ) / 3
      IRNH = LSEND
      LNXTRW = ( LIWK2 - LIRNH ) / 2
      NXTRW1 = IRNH + LIRNH
      NXTRW2 = NXTRW1 + LNXTRW
      IF ( NEWPT > LNXTRW .OR. LIRNH <= 0 ) GO TO 900

!  NXTROW( 1, . ) gives the link list. The list for column J starts
!                 in NXTROW( 1, J ) and ends when NXTROW( 1, K ) = - 1.
!  NXTROW( 2, . ) gives the position in H of the current link.

!  Initialize the link list which points to the row numbers which
!  are used in the columns of the assembled Hessian

      DO 20 I = 1, NUMVAR
         IWK( NXTRW1 + I ) = - 1
   20 CONTINUE

! -------------------------------------------------------
!  Form the rank-one second order term for the IG-th group
! -------------------------------------------------------

      DO 200 IG = 1, NG
         IF ( LOGI( GXEQX + IG ) ) GO TO 200
         IG1 = IG + 1
         LISTVS = IWK( LSTAGV + IG )
         LISTVE = IWK( LSTAGV + IG1 ) - 1

!  Form the J-th column of the rank-one matrix

         DO 190 L = LISTVS, LISTVE
            J = IWK( LSVGRP + L )
            IF ( J == 0 ) GO TO 190

!  Find the entry in row I of this column

            DO 180 K = LISTVS, LISTVE
               I = IWK( LSVGRP + K )
               IF ( I == 0 .OR. I > J ) GO TO 180

!  Obtain the appropriate storage location in H for the new entry

               ISTART = J
  150          CONTINUE
               INEXT = IWK( NXTRW1 + ISTART )
               IF ( INEXT == - 1 ) THEN
                  IF ( NEWPT > LNXTRW ) GO TO 900

!  The (I,J)-th location is empty. Place the new entry in this location
!  and add another link to the list

                  NNZH = NNZH + 1
                  IF ( NNZH > LIRNH ) GO TO 900
                  IWK( IRNH + NNZH ) = I
                  IWK( NXTRW1 + ISTART ) = NEWPT
                  IWK( NXTRW2 + ISTART ) = NNZH
                  IWK( NXTRW1 + NEWPT ) = - 1
                  NEWPT = NEWPT + 1
               ELSE

!  Continue searching the linked list for an entry in row I, column J

                  IF ( IWK( IRNH + IWK( NXTRW2 + ISTART ) )/=I ) THEN
                     ISTART = INEXT
                     GO TO 150
                  END IF
               END IF
  180       CONTINUE
  190    CONTINUE
  200 CONTINUE

! --------------------------------------------------------
!  Add on the low rank first order terms for the I-th group
! --------------------------------------------------------

      DO 300 IG = 1, NG
         IG1 = IG + 1

!  See if the group has any nonlinear elements

         DO 290 IELL = IWK( ISTADG + IG ), IWK( ISTADG + IG1 ) - 1
            IEL = IWK( IELING + IELL )
            LISTVS = IWK( ISTAEV + IEL )
            LISTVE = IWK( ISTAEV + IEL + 1 ) - 1
            DO 250 L = LISTVS, LISTVE
               J = IWK( IELVAR + L )
               IF ( J /= 0 ) THEN

!  The IEL-th element has an internal representation.
!  Compute the J-th column of the element Hessian matrix

!  Find the entry in row I of this column

                  DO 240 K = LISTVS, L
                     I = IWK( IELVAR + K )
                     IF ( I /= 0 ) THEN

!  Only the upper triangle of the matrix is stored

                        IF ( I <= J ) THEN
                           II = I
                           JJ = J
                        ELSE
                           II = J
                           JJ = I
                        END IF

!  Obtain the appropriate storage location in H for the new entry

                        ISTART = JJ
  230                   CONTINUE
                        INEXT = IWK( NXTRW1 + ISTART )
                        IF ( INEXT == - 1 ) THEN
                           IF ( NEWPT > LNXTRW ) GO TO 900

!  The (I,J)-th location is empty. Place the new entry in this location
!  and add another link to the list

                           NNZH = NNZH + 1
                           IF ( NNZH > LIRNH ) GO TO 900
                           IWK( IRNH + NNZH ) = II
                           IWK( NXTRW1 + ISTART ) = NEWPT
                           IWK( NXTRW2 + ISTART ) = NNZH
                           IWK( NXTRW1 + NEWPT ) = - 1
                           NEWPT = NEWPT + 1
                        ELSE

!  Continue searching the linked list for an entry in row I, column J

                           IF ( IWK( IRNH + IWK( NXTRW2 + ISTART ) ) &
 == II ) THEN
                           ELSE
                              ISTART = INEXT
                              GO TO 230
                           END IF
                        END IF
                     END IF
  240             CONTINUE
               END IF
  250       CONTINUE
  290    CONTINUE
  300 CONTINUE
      RETURN

!  Unsuccessful returns.

  900 CONTINUE
      WRITE( IOUT, 2000 )
      STOP

! Non-executable statements.

 2000 FORMAT( ' ** SUBROUTINE CDIMSH: Increase the size of IWK ' )

!  end of CDIMSH.

      END
