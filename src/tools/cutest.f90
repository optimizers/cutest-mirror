! THIS VERSION: CUTEST 1.0 - 24/10/2012 AT 13:00 GMT.

!-*-*-*-*-*-*-*-*-*-*-*-*-*- C U T E S T   M O D U l E -*-*-*-*-*-*-*-*-*-*-*-*-

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
      PUBLIC :: CUTEST_initialize_workspace, CUTEST_form_gradients,            &
                CUTEST_assemble_hessian, CUTEST_hessian_times_vector

!--------------------
!   P r e c i s i o n
!--------------------

      INTEGER, PARAMETER :: sp = KIND( 1.0E+0 )
      INTEGER, PARAMETER :: dp = KIND( 1.0D+0 )
      INTEGER, PARAMETER :: wp = dp

!----------------------
!   P a r a m e t e r s
!----------------------

!     INTEGER, PARAMETER :: lmin = 1
      INTEGER, PARAMETER :: lmin = 10000
      INTEGER, PARAMETER :: io_buffer = 75
      REAL ( KIND = wp ), PARAMETER :: zero = 0.0_wp
      REAL ( KIND = wp ), PARAMETER :: one = 1.0_wp

!-------------------------------------
!   G e n e r i c  i n t e r f a c e s
!-------------------------------------

!  define generic interfaces to routines for extending allocatable arrays

      INTERFACE CUTEST_extend_arrays
        MODULE PROCEDURE CUTEST_extend_array_real, CUTEST_extend_array_integer
      END INTERFACE

!-------------------------------------------------
!  D e r i v e d   t y p e   d e f i n i t i o n s
!-------------------------------------------------

!  =====================================
!  The CUTEST_assemble_type derived type
!  =====================================

      TYPE, PUBLIC :: CUTEST_assemble_type
        LOGICAL :: ptr_status
        INTEGER, DIMENSION( 30 ) :: ICNTL
        INTEGER, DIMENSION( 20 ) :: INFO
        REAL ( KIND = wp ), DIMENSION( 5 ) :: CNTL
      END TYPE CUTEST_assemble_type

!  =================================
!  The CUTEST_data_type derived type
!  =================================

      TYPE, PUBLIC :: CUTEST_data_type

!  Integer variables from the GLOBAL common block.

        INTEGER :: n, ng, nel, ntotel, nvrels, nnza, ngpvlu, nepvlu
        INTEGER :: ng1, nel1, nvargp, nvar2, nnonnz, nbprod, ntotin
 
        INTEGER :: lo, ch, liwork, lwork, ngng, la, lb, nobjgr, lu, lelvar
        INTEGER :: lstaev, lstadh, lntvar, lcalcf, leling, lintre, lft
        INTEGER :: lgxeqx, licna, lstada, lkndof, lgpvlu, lepvlu
        INTEGER :: lstadg, lgvals, lgscal, lescal, lvscal, lcalcg            
        INTEGER :: llink, lpos, lwtran, litran, l_link_e_u_v

        INTEGER :: lirnh = lmin
        INTEGER :: ljcnh = lmin
        INTEGER :: llink_min = lmin
        INTEGER :: lirnh_min = lmin
        INTEGER :: ljcnh_min = lmin
        INTEGER :: lh_min = lmin
        INTEGER :: lwtran_min = lmin
        INTEGER :: litran_min = lmin
        INTEGER :: lh = lmin
        INTEGER :: io_buffer = io_buffer

        LOGICAL :: alllin, densep, skipg
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
        LOGICAL, ALLOCATABLE, DIMENSION( : ) :: intrep
        LOGICAL, ALLOCATABLE, DIMENSION( : ) :: gxeqx
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

        INTEGER, ALLOCATABLE, DIMENSION( : ) :: ITRANS
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: LINK_col
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: POS_in_H
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: LINK_elem_uses_var
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: WTRANS

        INTEGER, ALLOCATABLE, DIMENSION( : ) :: ISYMMD
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: ISWKSP
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: ISTAJC
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: ISTAGV
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: ISVGRP
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: ISLGRP
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: IGCOLJ
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: IVALJR
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: IUSED 
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: ITYPER
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: ISSWTR
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: ISSITR
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: ISET
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: ISVSET
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: INVSET
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: IFREE
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: INDEX
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: IFREEC
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: INNONZ
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: LIST_elements
        INTEGER, ALLOCATABLE, DIMENSION( : , : ) :: ISYMMH
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: FUVALS_temp
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: P
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: X0
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: XCP
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: GX0
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: RADII
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DELTAX
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: QGRAD
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: GRJAC
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: GV_old
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : , : ) :: BND
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : , : ) :: BND_radius

        INTEGER, ALLOCATABLE, DIMENSION( : ) :: IW_asmbl
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: NZ_comp_w
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: W_ws
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: W_el
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: W_in
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: H_el
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: H_in

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

!  module procedures

    CONTAINS

!-  C U T E S T _ i n i t i a l i z e _ w o r k s p a c e  S U B R O U T I N E -

      SUBROUTINE CUTEST_initialize_workspace(                                  &
                        n, ng, nel, ntotel, nvrels, nnza, numvar, nvargp,      &
                        IELING, ISTADG, IELVAR, ISTAEV, INTVAR, ISTADH, ICNA,  &
                        ISTADA, ITYPEE, GXEQX, INTREP, altriv, direct, fdgrad, &
                        lfxi, lgxi, lhxi, lggfx, ldx, lnguvl, lnhuvl,  ntotin, &
                        ntype, nsets, maxsel, RANGE, iprint, iout, buffer,     &
! workspace
                        lwtran, litran, lwtran_min, litran_min, l_link_e_u_v,  &
                        llink_min, ITRANS, LINK_elem_uses_var, WTRANS,         &
                        ISYMMD, ISWKSP, ISTAJC, ISTAGV, ISVGRP, ISLGRP,        &
                        IGCOLJ, IVALJR, IUSED, ITYPER, ISSWTR, ISSITR,         &
                        ISET, ISVSET, INVSET, LIST_elements, ISYMMH,           &
                        IW_asmbl, NZ_comp_w, W_ws, W_el, W_in, H_el, H_in,     &
                        status, alloc_status, bad_alloc, skipg, KNDOFG )

!  Compute the starting addresses for the partitions of the workspace array
!  FUVALS. Also fill relevant portions of the workspace arrays WTRANS and ITRANS

!  History -
!   fortran 77 version originally released in CUTE, 20th June 1990
!   fortran 90 version released pre GALAHAD Version 1.0. February 1st 1995 as
!     INITW_initialize_workspace as part of the GALAHAD module INITW
!   fortran 2003 version released in CUTEst, 5th November 2012

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER, INTENT( IN ) :: n, ng, nel, ntotel, nvrels, nnza, numvar
      INTEGER, INTENT( IN ) :: iprint, iout, buffer
      INTEGER, INTENT( OUT ) :: lfxi, lgxi, lhxi, lggfx, ldx  
      INTEGER, INTENT( OUT ) :: lnguvl, lnhuvl, nvargp, status, alloc_status
      INTEGER, INTENT( OUT ) :: ntotin, ntype, nsets, maxsel
      LOGICAL, INTENT( IN ) :: direct, fdgrad, skipg
      LOGICAL, INTENT( OUT ) :: altriv
      CHARACTER ( LEN = 24 ), INTENT( OUT ) :: bad_alloc
      INTEGER, INTENT( IN ), DIMENSION( ntotel  ) :: IELING
      INTEGER, INTENT( IN ), DIMENSION( ng  + 1 ) :: ISTADA, ISTADG
      INTEGER, INTENT( IN ), DIMENSION( nel + 1 ) :: ISTAEV
      INTEGER, INTENT( IN ), DIMENSION( nvrels  ) :: IELVAR
      INTEGER, INTENT( IN ), DIMENSION( nnza    ) :: ICNA
      INTEGER, INTENT( OUT ), DIMENSION( nel + 1 ) :: ISTADH
      INTEGER, INTENT( INOUT ), DIMENSION( nel + 1 ) :: INTVAR
      INTEGER, INTENT( IN ), DIMENSION ( : ) :: ITYPEE
      LOGICAL, INTENT( IN ), DIMENSION( ng  ) :: GXEQX
      LOGICAL, INTENT( IN ), DIMENSION( nel ) :: INTREP

      INTEGER, INTENT( IN ), OPTIONAL, DIMENSION( ng ) :: KNDOFG

!-------------------------------------------------------------
!   D u m m y   A r g u m e n t s  f o r   w o r k s p a c e 
!-------------------------------------------------------------

      INTEGER, INTENT( INOUT ) :: lwtran, litran, lwtran_min, litran_min
      INTEGER, INTENT( INOUT ) :: l_link_e_u_v, llink_min

      INTEGER, ALLOCATABLE, DIMENSION( : ) :: ITRANS
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: LINK_elem_uses_var
      REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: WTRANS
  
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: ISYMMD
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: ISWKSP
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: ISTAJC
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: ISTAGV
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: ISVGRP
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: ISLGRP
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: IGCOLJ
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: IVALJR
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: IUSED 
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: ITYPER
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: ISSWTR
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: ISSITR
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: ISET
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: ISVSET
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: INVSET
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: LIST_elements
      INTEGER, ALLOCATABLE, DIMENSION( : , : ) :: ISYMMH
  
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: IW_asmbl
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: NZ_comp_w
      REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: W_ws
      REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: W_el
      REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: W_in
      REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: H_el
      REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: H_in

!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------

      INTERFACE
        SUBROUTINE RANGE ( ielemn, transp, W1, W2, nelvar, ninvar, ieltyp,     &
                           lw1, lw2 )
        INTEGER, INTENT( IN ) :: ielemn, nelvar, ninvar, ieltyp, lw1, lw2
        LOGICAL, INTENT( IN ) :: transp
        REAL ( KIND = KIND( 1.0D+0 ) ), INTENT( IN ), DIMENSION ( lw1 ) :: W1
        REAL ( KIND = KIND( 1.0D+0 ) ), INTENT( OUT ), DIMENSION ( lw2 ) :: W2
        END SUBROUTINE RANGE
      END INTERFACE

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: i, j, k, l, iielts, ientry, ig, is
      INTEGER :: nsizeh, nel1, iel, lwfree, lifree, lnwksp
      INTEGER :: nelvr, liwfro, lwfreo, iell, itype, isofar
      INTEGER :: istarj, ivarp1, ivar, jset, inext, newvar, newset
      INTEGER :: ipt, istrt, ninvr, ii, kk, ll
      INTEGER :: nwtran, mwtran, uwtran, llink, mlink, nlink, ulink
      INTEGER :: nitran, mitran, uitran, maxsin
      LOGICAL :: alllin, vrused, reallocate

!     CHARACTER ( LEN = 80 ) :: array

!  Set constants

      nel1 = nel + 1
      alllin = nel == 0

!  Set up INTVAR, the starting addresses for the element gradients with
!  respect to their internal variables. Also compute maxsin, the maximum
!  number of internal variables in an element

      IF ( .NOT. alllin ) THEN
        k = INTVAR( 1 )
        maxsin = k
        INTVAR( 1 ) = nel1
        DO iel = 2, nel
          l = INTVAR( iel )
          INTVAR( iel ) = INTVAR( iel - 1 ) + k
          k = l
          maxsin = MAX( maxsin, k )
        END DO
        INTVAR( nel1 ) = INTVAR( nel ) + k
      ELSE
        INTVAR( 1 ) = 1
        maxsin = 0
      END IF

!  Compute the total number of internal variables

      ntotin = INTVAR( nel1 ) - INTVAR( 1 )

!  Calculate the length, iielts, of workspace required to determine which
!  elements use each of the variables. Also find the maximum number of
!  variables in an element, maxsel. This is a dummy run merely to calculate
!  the space required

      llink = n + llink_min
      reallocate = .TRUE.
      IF ( ALLOCATED( LINK_elem_uses_var ) ) THEN
        IF ( SIZE( LINK_elem_uses_var ) < llink ) THEN
          DEALLOCATE( LINK_elem_uses_var )
        ELSE ; llink = SIZE( LINK_elem_uses_var ) ; reallocate = .FALSE.
        END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( LINK_elem_uses_var( llink ), STAT = alloc_status )
        IF ( alloc_status /= 0 ) THEN 
          bad_alloc = 'LINK_e' ; GO TO 600 ; END IF
      END IF

!  LINK_elem_uses_var( i ) will be used as a list of links chaining the 
!  elements using variable i. If LINK_elem_uses_var( i ) is negative, the 
!  list is empty

      LINK_elem_uses_var( : n ) = - 1
      iielts = n ; maxsel = 0
      IF ( .NOT. alllin ) THEN

!  Loop over the groups, considering each nonlinear element in turn

        DO iel = 1, nel
          maxsel = MAX( maxsel, ISTAEV( iel + 1 ) - ISTAEV( iel ) )
        END DO
        DO i = 1, ntotel
          iel = IELING( i )

!  Loop on the variables from the I-th element

          DO k = ISTAEV( iel ), ISTAEV( iel + 1 ) - 1
            ientry = IELVAR( k )
            IF ( LINK_elem_uses_var( ientry ) >= 0 ) THEN

!  If we have reached the end of the list of the elements using the variable
!  IELVAR( K ), add the IEL-th element to it. Otherwise, find the next entry
!  in the list

  30          CONTINUE
              IF ( LINK_elem_uses_var( ientry ) > 0 ) THEN
                ientry = LINK_elem_uses_var( ientry )
                GO TO 30
              ELSE
               IF ( iielts == llink ) THEN
                  nlink = llink
                  ulink = iielts; mlink = iielts + 1
                  CALL CUTEST_extend_array( LINK_elem_uses_var, llink, ulink,  &
                                            nlink, mlink, buffer, status,      &
                                            alloc_status)
                  IF ( status /= 0 ) THEN
                    bad_alloc = 'LINK_elem_uses_var' ; GO TO 610 ; END IF
                  llink = nlink
                END IF
                iielts = iielts + 1
                LINK_elem_uses_var( ientry ) = iielts
                LINK_elem_uses_var( iielts ) = 0
              END IF
            ELSE

!  The list of elements involving the variable IELVAR( K ) was
!  previously empty. Indicate that the list has now been started and
!  that its end has been reached

              LINK_elem_uses_var( ientry ) = 0
            END IF
          END DO
        END DO
      END IF
       
      l_link_e_u_v = iielts

!  -- Calculate the starting addresses for the integer workspace --

!  ISWKSP( j ), j = 1, ..., MAX( ntotel, nel, n + n ), is used for
!  workspace by the matrix-vector product subroutine HSPRD

      IF ( direct ) THEN
        lnwksp = MAX( MAX( ntotel, nel ), n + n )
      ELSE
        lnwksp = MAX( MAX( ntotel, nel ), n )
      END IF
      
      reallocate = .TRUE.
      IF ( ALLOCATED( ISWKSP ) ) THEN
        IF ( SIZE( ISWKSP ) < lnwksp ) THEN ; DEALLOCATE( ISWKSP )
        ELSE ; reallocate = .FALSE. ; END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( ISWKSP( lnwksp ), STAT = alloc_status )
        IF ( alloc_status /= 0 ) THEN 
          bad_alloc = 'ISWKSP' ; GO TO 600 ; END IF
      END IF

!  IUSED( j ), j = 1, ..., MAX( n, ng ) Will be used as workspace by
!  the matrix-vector product subroutine HSPRD

      reallocate = .TRUE.
      IF ( ALLOCATED( IUSED ) ) THEN
        IF ( SIZE( IUSED ) < MAX( n, ng ) ) THEN ; DEALLOCATE( IUSED )
        ELSE ; reallocate = .FALSE. ; END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( IUSED( MAX( n, ng ) ), STAT = alloc_status )
        IF ( alloc_status /= 0 ) THEN 
          bad_alloc = 'IUSED' ; GO TO 600 ; END IF
      END IF

!  ISLGRP( j ), j = 1, ..., ntotel, will contain the number of the group
!  which uses nonlinear element j

      reallocate = .TRUE.
      IF ( ALLOCATED( ISLGRP ) ) THEN
        IF ( SIZE( ISLGRP ) < ntotel ) THEN ; DEALLOCATE( ISLGRP )
        ELSE ; reallocate = .FALSE. ; END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( ISLGRP( ntotel ), STAT = alloc_status )
        IF ( alloc_status /= 0 ) THEN 
          bad_alloc = 'ISLGRP' ; GO TO 600 ; END IF
      END IF

!  ISTAJC( j ), j = 1, ..., n, will contain the starting addresses for
!  the list of nontrivial groups which use the j-th variable.
!  ISTAJC( n + 1 ) will point to the first free location in IGCOLJ after
!  the list of nontrivial groups for the n-th variable

      reallocate = .TRUE.
      IF ( ALLOCATED( ISTAJC ) ) THEN
        IF ( SIZE( ISTAJC ) < n + 1 ) THEN ; DEALLOCATE( ISTAJC )
        ELSE ; reallocate = .FALSE. ; END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( ISTAJC( n + 1 ), STAT = alloc_status )
        IF ( alloc_status /= 0 ) THEN 
          bad_alloc = 'ISTAJC' ; GO TO 600 ; END IF
      END IF

!  ISTAGV( j ), j = 1, ..., ng, will contain the starting addresses for
!  the list of variables which occur in the J-th group. ISTAGV( ng + 1 )
!  will point to the first free location in ISVGRP after the list of variables
!  for the NG-th group

      reallocate = .TRUE.
      IF ( ALLOCATED( ISTAGV ) ) THEN
        IF ( SIZE( ISTAGV ) < ng + 1 ) THEN ; DEALLOCATE( ISTAGV )
        ELSE ; reallocate = .FALSE. ; END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( ISTAGV( ng + 1 ), STAT = alloc_status )
        IF ( alloc_status /= 0 ) THEN 
          bad_alloc = 'ISTAGV' ; GO TO 600 ; END IF
      END IF

!  Allocate LIST_elements

      reallocate = .TRUE.
      IF ( ALLOCATED( LIST_elements ) ) THEN
        IF ( SIZE( LIST_elements ) < l_link_e_u_v ) THEN 
          DEALLOCATE( LIST_elements )
        ELSE ; reallocate = .FALSE. ; END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( LIST_elements( l_link_e_u_v ), STAT = alloc_status )
        IF ( alloc_status /= 0 ) THEN 
          bad_alloc = 'LIST_e' ; GO TO 600 ; END IF
      END IF

!  Determine which elements use each variable. Initialization

      IF ( .NOT. alllin ) THEN

!  LINK_elem_uses_var( i ) will be used as a list of links chaining the 
!  elements using variable i. If LINK_elem_uses_var( i ) is negative, the 
!  list is empty

        LINK_elem_uses_var( : n ) = - 1
        LIST_elements( : n ) = - 1   ! needed for epcf90 debugging compiler
        iielts = n

!  Loop over the groups, considering each nonlinear element in turn

        DO i = 1, ntotel
          iel = IELING( i )

!  Loop on the variables of the I-th element

          DO k = ISTAEV( iel ), ISTAEV( iel + 1 ) - 1
            ientry = IELVAR( k )
            IF ( LINK_elem_uses_var( ientry ) >= 0 ) THEN

!  If we have reached the end of the list of the elements using the variable
!  IELVAR( K ), add the I-th element to it and record that the end of the list
!  has occured. Otherwise, find the next entry in the list

  110         CONTINUE
              IF ( LINK_elem_uses_var( ientry ) > 0 ) THEN
                ientry = LINK_elem_uses_var( ientry )
                GO TO 110
              ELSE
                iielts = iielts + 1
                LINK_elem_uses_var( ientry ) = iielts
                LINK_elem_uses_var( iielts ) = 0
                LIST_elements( iielts ) = i
              END IF
            ELSE

!  The list of elements involving the variable IELVAR( K ) was previously
!  empty. Indicate that the list has now been started, record the element
!  which contains the variable and indicate that the end of the list has
!  been reached

              LINK_elem_uses_var( ientry ) = 0
              LIST_elements( ientry ) = i
            END IF
          END DO
        END DO
      END IF

!  Set up symmetric addresses for the upper triangular storage
!  schemes for the element hessians

      IF ( maxsin > 0 ) THEN
        reallocate = .TRUE.
        IF ( ALLOCATED( ISYMMD ) ) THEN
           IF ( SIZE( ISYMMD ) < maxsin ) THEN ; DEALLOCATE( ISYMMD )
           ELSE ; reallocate = .FALSE.
           END IF
        END IF
        IF ( reallocate ) THEN 
           ALLOCATE( ISYMMD( maxsin ), STAT = alloc_status )
           IF ( alloc_status /= 0 ) THEN ; bad_alloc = 'ISYMMD' ; GO TO 600
           END IF
        END IF
        
        reallocate = .TRUE.
        IF ( ALLOCATED( ISYMMH ) ) THEN
          IF ( SIZE( ISYMMH, 1 ) /= maxsin .OR. SIZE( ISYMMH, 2 ) /= maxsin ) &
            THEN  ; DEALLOCATE( ISYMMH ) ; ELSE ; reallocate = .FALSE. ; END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( ISYMMH( maxsin, maxsin ), STAT = alloc_status )
          IF ( alloc_status /= 0 ) THEN ; bad_alloc = 'ISYMMH' ; GO TO 600
          END IF
        END IF
        
        CALL CUTEST_symmh( maxsin, ISYMMH, ISYMMD )
      ELSE
        ALLOCATE( ISYMMD( 0 ), STAT = alloc_status )
        IF ( alloc_status /= 0 ) THEN ; bad_alloc = 'ISYMMD' ; GO TO 600
        END IF
        ALLOCATE( ISYMMH( 0, 0 ), STAT = alloc_status )
        IF ( alloc_status /= 0 ) THEN ; bad_alloc = 'ISYMMH' ; GO TO 600
        END IF
      END IF

!  Set up the starting addresses for the element Hessians
!  with respect to their internal variables and a pointer beyond
!  the end of the space required for the Hessians

      lggfx = INTVAR( nel1 )
      IF ( .NOT. alllin ) THEN
        DO i = 1, nel
          ISTADH( i ) = lggfx
          nsizeh = INTVAR( i + 1 ) - INTVAR( i )
          lggfx = lggfx + nsizeh * ( nsizeh + 1 ) / 2
        END DO
      END IF
      ISTADH( nel1 ) = lggfx

!  ALTRIV specifies whether all the groups are trivial

      altriv = .TRUE.

!  Pass 1: Count the total number of variables in all the groups, nvargp

      nvargp = 0

!  Start by initializing the counting array to zero

      ISWKSP( : numvar ) = 0

!  Loop over the groups. See if the IG-th group is trivial

      DO ig = 1, ng

!  Check to see if all of the groups are trivial

        IF ( skipg ) THEN ; IF ( KNDOFG( ig ) == 0 ) CYCLE ; END IF
        IF ( .NOT. GXEQX( ig ) ) altriv = .FALSE.

!  Loop over the nonlinear elements from the IG-th group

        DO k = ISTADG( ig ), ISTADG( ig + 1 ) - 1
          iel = IELING( k )

!  Run through all the elemental variables changing the I-th entry of
!  ISWKSP from zero to one if variable I appears in an element

          DO j = ISTAEV( iel ), ISTAEV( iel + 1 ) - 1
            i = IELVAR( j )
            IF ( ISWKSP( i ) < ig ) THEN
              ISWKSP( i ) = ig
              nvargp = nvargp + 1
            END IF
          END DO
        END DO

!  Consider variables which arise from the linear element

        DO j = ISTADA( ig ), ISTADA( ig + 1 ) - 1
          i = ICNA( j )
          IF ( i <= numvar ) THEN
            IF ( ISWKSP( i ) < ig ) THEN
               ISWKSP( i ) = ig
               nvargp = nvargp + 1
            END IF
          END IF
        END DO
      END DO

!  ISVGRP( j ), j = 1, ..., nvargp, will contain the indices of the
!  variables which are used by each group in turn. Those for group i occur
!  in locations ISTAGV( i ) to ISTAGV( i + 1 ) - 1

!  Allocate the array ISVGRP

      reallocate = .TRUE.
      IF ( ALLOCATED( ISVGRP ) ) THEN
        IF ( SIZE( ISVGRP ) < nvargp ) THEN ; DEALLOCATE( ISVGRP )
        ELSE ; reallocate = .FALSE. ; END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( ISVGRP( nvargp ), STAT = alloc_status )
        IF ( alloc_status /= 0 ) THEN 
          bad_alloc = 'ISVGRP' ; GO TO 600 ; END IF
      END IF

!  Store the indices of variables which appears in each group and how many
!  groups use each variable. Reinitialize counting arrays to zero

      ISTAJC( 2 : n + 1 ) = 0
      ISWKSP( : numvar ) = 0

!  Pass 2: store the list of variables

      nvargp = 0
      ISTAGV( 1 ) = 1

!  Loop over the groups. See if the IG-th group is trivial

      DO ig = 1, ng

        IF ( skipg ) THEN 
          IF ( KNDOFG( ig ) == 0 ) THEN
            ISLGRP( ISTADG( ig ) : ISTADG( ig + 1 ) - 1 ) = ig
            ISTAGV( ig + 1 ) = nvargp + 1
            CYCLE
          END IF
        END IF

!  Again, loop over the nonlinear elements from the IG-th group

        DO k = ISTADG( ig ), ISTADG( ig + 1 ) - 1
          iel = IELING( k )

!  Run through all the elemental variables changing the I-th entry of
!  ISWKSP from zero to one if variable I appears in an element

          DO j = ISTAEV( iel ), ISTAEV( iel + 1 ) - 1
            i = IELVAR( j )
            IF ( ISWKSP( i ) < ig ) THEN
              ISWKSP( i ) = ig

!  Record the nonlinear variables from the ig-th group

              nvargp = nvargp + 1
              ISVGRP( nvargp ) = i
            END IF
          END DO

!  Record that nonlinear element K occurs in group IELGRP( IEL )

          ISLGRP( k ) = ig
        END DO

!  Consider variables which arise from the linear element

        DO j = ISTADA( ig ), ISTADA( ig + 1 ) - 1
          i = ICNA( j )
          IF ( i <= numvar ) THEN
            IF ( ISWKSP( i ) < ig ) THEN
              ISWKSP( i ) = ig

!  Record the linear variables from the ig-th group

              nvargp = nvargp + 1
              ISVGRP( nvargp ) = i
            END IF
          END IF
        END DO

!  Record that one further nontrivial group uses variable l - 1

        IF ( .NOT. GXEQX( ig ) ) THEN
          DO j = ISTAGV( ig ), nvargp
            l = ISVGRP( j ) + 1
            ISTAJC( l ) = ISTAJC( l ) + 1
          END DO
        END IF

!  Record the starting address of the variables in the next group

        ISTAGV( ig + 1 ) = nvargp + 1
      END DO
      ISWKSP( : n ) = 0

!  IGCOLJ( j ), j = 1, ..., nvargp, will contain the indices of the
!  nontrivial groups which use each variable in turn. Those for variable i
!  occur in locations ISTAJC( i ) to ISTAJC( i + 1 ) - 1

      reallocate = .TRUE.
      IF ( ALLOCATED( IGCOLJ ) ) THEN
        IF ( SIZE( IGCOLJ ) < nvargp ) THEN ; DEALLOCATE( IGCOLJ )
        ELSE ; reallocate = .FALSE. ; END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( IGCOLJ( nvargp ), STAT = alloc_status )
        IF ( alloc_status /= 0 ) THEN 
          bad_alloc = 'IGCOLJ' ; GO TO 600 ; END IF
      END IF

!  IVALJR( j ), j = 1, ..., nvargp, will contain the positions in GRJAC
!  of the nonzeros of the Jacobian
!  of the groups corresponding to the variables as ordered in ISVGRP( j )

      reallocate = .TRUE.
      IF ( ALLOCATED( IVALJR ) ) THEN
        IF ( SIZE( IVALJR ) < nvargp ) THEN ; DEALLOCATE( IVALJR )
        ELSE ; reallocate = .FALSE. ; END IF
      END IF
      IF ( reallocate ) THEN 
        ALLOCATE( IVALJR( nvargp ), STAT = alloc_status )
        IF ( alloc_status /= 0 ) THEN ; bad_alloc = 'IVALJR' ; GO TO 600
        END IF
      END IF

!  Set the starting addresses for the lists of nontrivial groups which use
!  each variable in turn

      k = 1
      ISTAJC( k ) = 1
      DO i = 2, n + 1
        k = k + 1
        ISTAJC( k ) = ISTAJC( k ) + ISTAJC( k - 1 )
      END DO

!  Consider the IG-th group in order to associate variables with groups

      DO ig = 1, ng
        IF ( skipg ) THEN ; IF ( KNDOFG( ig ) == 0 ) CYCLE ; END IF
        IF ( .NOT. GXEQX( ig ) ) THEN
          DO i = ISTAGV( ig ), ISTAGV( ig + 1 ) - 1
            l = ISVGRP( i )

!  Record that group IG uses variable ISVGRP( i )

            j = ISTAJC( l )
            IGCOLJ( j ) = ig

!  Store the locations in the Jacobian of the groups of the nonzeros
!  corresponding to each variable in the IG-TH group. Increment the starting
!  address for the pointer to the next group using variable ISVGRP( i )

            IVALJR( i ) = j
            ISTAJC( l ) = j + 1
          END DO
        END IF
      END DO

!  Reset the starting addresses for the lists of groups using each variable

      DO i = n, 2, - 1
        ISTAJC( i ) = ISTAJC( i - 1 )
      END DO
      ISTAJC( 1 ) = 1

!  Initialize workspace values for subroutine HSPRD

      IUSED( : MAX( n, ng ) ) = 0

!  Initialize general workspace arrays

      maxsin = MAX( 1, maxsin )
      maxsel = MAX( 1, maxsel )

      IF ( ALLOCATED( IW_asmbl ) ) DEALLOCATE( IW_asmbl )
      ALLOCATE( IW_asmbl( n ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN ; bad_alloc = 'IW_asmb' ; GO TO 600
      END IF
      
      IF ( ALLOCATED( NZ_comp_w ) ) DEALLOCATE( NZ_comp_w )
      ALLOCATE( NZ_comp_w( ng ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN ; bad_alloc = 'NZ_com' ; GO TO 600
      END IF
      
      IF ( ALLOCATED( W_ws ) ) DEALLOCATE( W_ws )
      ALLOCATE( W_ws( MAX( n, ng ) ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN ; bad_alloc = 'W_ws' ; GO TO 600
      END IF
      
      IF ( ALLOCATED( W_el ) ) DEALLOCATE( W_el )
      ALLOCATE( W_el( maxsel ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN ; bad_alloc = 'W_el' ; GO TO 600
      END IF
      
      IF ( ALLOCATED( W_in ) ) DEALLOCATE( W_in )
      ALLOCATE( W_in( maxsin ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN ; bad_alloc = 'W_in' ; GO TO 600
      END IF
      
      IF ( ALLOCATED( H_el ) ) DEALLOCATE( H_el )
      ALLOCATE( H_el( maxsel ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN ; bad_alloc = 'H_el' ; GO TO 600
      END IF
      
      IF ( ALLOCATED( H_in ) ) DEALLOCATE( H_in )
      ALLOCATE( H_in( maxsin ), STAT = alloc_status )
      IF ( alloc_status /= 0 ) THEN ; bad_alloc = 'H_in' ; GO TO 600
      END IF

!  Define further partitions of the workspace whenever finite-difference
!  gradients are used

      IF ( fdgrad ) THEN

!  The range transformation for each nonlinear element is of a given type.
!  Suppose there are NTYPE non-trivial types. ITYPER( i ) gives the type
!  of nonlinear element i for i = 1, ...., nel

        reallocate = .TRUE.
        IF ( ALLOCATED( ITYPER ) ) THEN
          IF ( SIZE( ITYPER ) < nel ) THEN ; DEALLOCATE( ITYPER )
          ELSE ; reallocate = .FALSE. ; END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( ITYPER( nel ), STAT = alloc_status )
          IF ( alloc_status /= 0 ) THEN
            bad_alloc = 'ITYPER' ; GO TO 600 ; END IF
        END IF

!  The range transformation from elemental to internal variables is defined by
!  a matrix W. For each non-trivial transformation, the matrix W is recorded.
!  The information for the I-th type starts in location 
!  ISSWTR( i ), i = 1, ...., ntype of WTRANS

        reallocate = .TRUE.
        IF ( ALLOCATED( ISSWTR ) ) THEN
          IF ( SIZE( ISSWTR ) < nel ) THEN ; DEALLOCATE( ISSWTR )
          ELSE ; reallocate = .FALSE. ; END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( ISSWTR( nel ), STAT = alloc_status )
          IF ( alloc_status /= 0 ) THEN
            bad_alloc = 'ISSWTR' ; GO TO 600 ; END IF
        END IF

!  For each type of nonlinear element using a nontrivial range transformation,
!  integer information is also recorded. The information for the i-th type
!  starts in location ISSITR( i ), i = 1, ...., ntype of ITRANS

        reallocate = .TRUE.
        IF ( ALLOCATED( ISSITR ) ) THEN
          IF ( SIZE( ISSITR ) < nel ) THEN ; DEALLOCATE( ISSITR )
          ELSE ; reallocate = .FALSE. ; END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( ISSITR( nel ), STAT = alloc_status )
          IF ( alloc_status /= 0 ) THEN 
            bad_alloc = 'ISSITR' ; GO TO 600 ; END IF
        END IF

!  The following pieces of integer information are recorded for the I-th type
!  of nonlinear element:

!    ITRANS( ISSITR( i ) + 1 ):
!            the number of internal variables, ninvr
!    ITRANS( ISSITR( i ) + 2 ):
!            the number of elemental variables, nelvr
!    ITRANS( ISSITR( i ) + 2 + j ), j = 1, ..., ninvr + nelvr:
!            pivot sequences for the LU factors of W.

!  After the factorization and compression, only ninvr linearly independent
!  columns of W are stored

!  -- Make an initial allocation of WTRANS and ITRANS

        reallocate = .TRUE.
        IF ( ALLOCATED( WTRANS ) ) THEN
          IF ( SIZE( WTRANS ) < lwtran_min ) THEN
             DEALLOCATE( WTRANS ) ; lwtran = lwtran_min
          ELSE ; reallocate = .FALSE. ; END IF
        ELSE ; lwtran = lwtran_min ; END IF
        IF ( reallocate ) THEN 
          ALLOCATE( WTRANS( lwtran ), STAT = alloc_status )
          IF ( alloc_status /= 0 ) THEN 
             bad_alloc = 'WTRANS' ; GO TO 600 ; END IF
        END IF
         
        reallocate = .TRUE.
        IF ( ALLOCATED( ITRANS ) ) THEN
          IF ( SIZE( ITRANS ) < litran_min ) THEN
            DEALLOCATE( ITRANS ) ; litran = litran_min
          ELSE ; reallocate = .FALSE. ; END IF
        ELSE ; litran = litran_min ; END IF
        IF ( reallocate ) THEN 
          ALLOCATE( ITRANS( litran ), STAT = alloc_status )
          IF ( alloc_status /= 0 ) THEN 
            bad_alloc = 'ITRANS' ; GO TO 600 ; END IF
        END IF

!  ---------------------------------------------------
!  Consider only elements which use internal variables
!  ---------------------------------------------------

        ntype = 0 ; lwfree = 1 ; lifree = 1 
        W_el( : maxsel ) = zero

!  Loop over all nonlinear elements

  LIEL: DO iel = 1, nel
          IF ( INTREP( iel ) ) THEN

!  Calculate the range transformation matrix WTRANS

            is = ISTAEV( iel )
            ninvr = INTVAR( iel + 1 ) - INTVAR( iel )
            nelvr = ISTAEV( iel + 1 ) - is
            mwtran = lwfree + ninvr * nelvr - 1
                  
!  Ensure that there is enough space

            IF ( mwtran > lwtran ) THEN
              nwtran = 2 * ( lwfree + ninvr * nelvr - 1 )
              uwtran = lwfree - 1
              CALL CUTEST_extend_array( WTRANS, lwtran, uwtran, nwtran,       &
                                        mwtran, buffer, status, alloc_status )
              IF ( status /= 0 ) THEN
                bad_alloc = 'WTRANS' ; GO TO 610 ; END IF
              lwtran = nwtran
            END IF
  
            k = lwfree - 1
            is = is - 1
  
            DO i = 1, nelvr
              W_el( i ) = one
              CALL RANGE ( iel, .FALSE., W_el( : nelvr ),                     &
                           WTRANS( k + 1 : k + ninvr ), nelvr, ninvr,         &
                           ITYPEE( iel ), nelvr, ninvr )
              W_el( i ) = zero
              k = k + ninvr

!  Check to see if any of the columns belong to duplicated variables

              ii = IELVAR( is + i )
              DO j = 1, i - 1
                IF ( IELVAR( is + j ) == ii ) GO TO 300
              END DO
              CYCLE

!  Amalgamate columns from duplicate variables

  300         CONTINUE
              kk = lwfree + ( j - 1 ) * ninvr - 1
              ll = k - ninvr
              WTRANS( kk + 1 : kk + ninvr ) =                                 &
                WTRANS( kk + 1 : kk + ninvr ) + WTRANS( ll + 1 : ll + ninvr )
              WTRANS( ll + 1 : ll + ninvr ) = zero
            END DO

!  Compare this transformation matrix with previous ones

        LI: DO i = 1, ntype
              IF ( ITRANS( ISSITR( i ) ) /= ninvr .OR.                        &
                   ITRANS( ISSITR( i ) + 1 ) /= nelvr ) CYCLE LI
              DO j = 0, ninvr * nelvr - 1
                IF ( WTRANS( lwfree + j ) /=                                  &
                     WTRANS( ISSWTR( i ) + j ) ) CYCLE LI
              END DO

!  The transformation is an existing one. Record which.

              ITYPER( iel ) = i
              CYCLE LIEL
            END DO LI

            mitran = lifree + ninvr + nelvr + 1
                  
!  Ensure that there is enough space

            IF ( mitran > litran ) THEN
              nitran = 2 * ( lifree + ninvr + nelvr + 1 )
              uitran = lifree - 1
              CALL CUTEST_extend_array( ITRANS, litran, uitran, nitran,       &
                                        mitran, buffer, status, alloc_status )
              IF ( status /= 0 ) THEN 
                bad_alloc = 'ITRANS' ; GO TO 610 ; END IF
              litran = nitran
            END IF

!  The transformation defines a new type. Record its details

            ntype = ntype + 1
            ITYPER( iel ) = ntype
            ITRANS( lifree ) = ninvr
            ITRANS( lifree + 1 ) = nelvr
            ITRANS( lifree + 2 : mitran + 1 ) = 0
            ISSITR( ntype ) = lifree
            ISSWTR( ntype ) = lwfree
            lifree = lifree + 2 + ninvr + nelvr
            lwfree = lwfree + ninvr * nelvr
          ELSE
            ITYPER( iel ) = 0
          END IF
        END DO LIEL

!  For each type of element with internal variables:

        DO i = 1, ntype
          lwfree = ISSWTR( i ) ; lifree = ISSITR( i )
          ninvr = ITRANS( lifree )
          nelvr = ITRANS( lifree + 1 )

!  Factorize W. Use Gaussian elimination with complete pivoting.
!  Determine the "most independent" set of columns of W

          CALL CUTEST_gauss_elim(                                             &
              ninvr, nelvr, ITRANS( lifree + 2 : lifree + ninvr + 1 ),        &
              ITRANS( lifree + ninvr + 2 : lifree + ninvr + nelvr + 1 ),      &
              WTRANS( lwfree : lwfree + ninvr * nelvr - 1 ) )
     
        END DO

!  Compress the data structures to remove redundant information

!  Compress integer data

        litran = 0
        lwtran = 0
        DO i = 1, ntype
          liwfro = ISSITR( i ) - 1
          ninvr = ITRANS( liwfro + 1 )
          k = 2 * ( ninvr + 1 )
          DO j = 1, k
             ITRANS( litran + j ) = ITRANS( liwfro + j )
          END DO
          ISSITR( i ) = litran + 1
          litran = litran + k

!  Compress real data

          lwfreo = ISSWTR( i ) - 1
          DO j = 1, ninvr * ninvr
             WTRANS( lwtran + j ) = WTRANS( lwfreo + j )
          END DO
          ISSWTR( i ) = lwtran + 1
          lwtran = lwtran + ninvr * ninvr
        END DO

!  ----------------------------------------------------------------------
!  The list of variables is allocated to nsets disjoints sets. Variable I
!  occurs in set ISET
!  ----------------------------------------------------------------------

        reallocate = .TRUE.
        IF ( ALLOCATED( ISET ) ) THEN
          IF ( SIZE( ISET ) < n ) THEN ; DEALLOCATE( ISET )
          ELSE ; reallocate = .FALSE. ; END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( ISET( n ), STAT = alloc_status )
          IF ( alloc_status /= 0 ) THEN 
            bad_alloc = 'ISET' ; GO TO 600 ; END IF
        END IF
        
        reallocate = .TRUE.
        IF ( ALLOCATED( ISVSET ) ) THEN
          IF ( SIZE( ISVSET ) < n + 2 ) THEN ; DEALLOCATE( ISVSET )
          ELSE ; reallocate = .FALSE. ; END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( ISVSET( n + 2 ), STAT = alloc_status )
          IF ( alloc_status /= 0 ) THEN 
            bad_alloc = 'ISVSET' ; GO TO 600 ; END IF
        END IF
        
        reallocate = .TRUE.
        IF ( ALLOCATED( INVSET ) ) THEN
          IF ( SIZE( INVSET ) < MAX( n + 1, nel ) ) THEN ; DEALLOCATE( INVSET )
          ELSE ; reallocate = .FALSE. ; END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( INVSET( MAX( n + 1, nel ) ), STAT = alloc_status )
          IF ( alloc_status /= 0 ) THEN 
            bad_alloc = 'INVSET' ; GO TO 600 ; END IF
        END IF

!  Assign initial set numbers to each variable

        nsets = 0
        ISET( : n ) = n

!  Use the Curtis-Powell-Reid algorithm to determine which set each variable
!  belongs to. Loop over the variables.

        DO i = 1, n

!  Loop over the elements which use variable i. The elements are obtained from
!  a linked-list

          vrused = .FALSE.
          ipt = LINK_elem_uses_var( i )
          IF ( ipt >= 0 ) THEN
            iell = LIST_elements( i )
  420       CONTINUE
            iel = IELING( iell )
            itype = ITYPER( iel )

!  Check that the variable belongs to the "independence" set of elements with
!  internal variables

            IF ( itype > 0 ) THEN
              lifree = ISSITR( itype )
              ninvr = ITRANS( lifree )
              DO j = 1, ninvr
                k = j - 1
                l = ITRANS( lifree + ninvr + 1 + j ) - 1
                IF ( i == IELVAR( ISTAEV( iel ) + l ) ) GO TO 440
              END DO
              GO TO 450
  440         CONTINUE
            END IF
            vrused = .TRUE.
  450       CONTINUE

!  Loop over the complete list of variables used by element iel

!DIR$ IVDEP
            DO j = ISTAEV( iel ), ISTAEV( iel + 1 ) - 1

!  If variable IV is used, flag the set which contains it

              INVSET( ISET( IELVAR( j ) ) ) = 1
            END DO

!  Check the link-list to see if further elements use the variable

            IF ( ipt > 0 ) THEN
              iell = LIST_elements( ipt )
              ipt = LINK_elem_uses_var( ipt )
              GO TO 420
            END IF
          END IF

!  See if the variable may be placed in the first nsets sets

          IF ( vrused ) THEN
            DO j = 1, nsets
              IF ( INVSET( j ) == 0 ) GO TO 480
              INVSET( j ) = 0
            END DO

!  The variable needs a new set

            nsets = nsets + 1

!  The variable will be placed in set j

            j = nsets

  480       CONTINUE
            ISET( i ) = j

!  Reset the flags to zero

            INVSET( j : nsets ) = 0
          ELSE

!  The variable is not to be used

            ISET( i ) = n
          END IF
        END DO

!  Check that there is at least one set

        IF ( nsets /= 0 ) THEN
          ISET( : n ) = MIN( ISET( : n ), nsets + 1 )

!  ---------------------------------------------------------------------
!  Obtain a list, INVSET, of the variables corresponding to each set
!  ---------------------------------------------------------------------

!  Clear ISVSET

          ISVSET( 2 : nsets + 2 ) = 0

!  Count the number of elements in each set and store in ISVSET.
!  Negate the set numbers in ISET, so that they are flagged as
!  ISET is gradually overwritten by variable indices.

          DO k = 1, n
            j = ISET( k )
            ISET( k ) = - j
            ISVSET( j + 1 ) = ISVSET( j + 1 ) + 1
          END DO

!  Compute the starting addresses for each set within IISET

          ISVSET( 1 ) = 1
          DO j = 2, nsets + 2
            ISVSET( j ) = ISVSET( j ) + ISVSET( j - 1 )
          END DO

!  Store in INVSET the variable whose set number is the
!  ISVSET( j )-th entry of INVSET

          isofar = 0
          DO j = 1, nsets + 1
            istarj = ISVSET( j )
            DO ivarp1 = isofar + 1, n
               IF ( istarj < ivarp1 ) GO TO 530
            END DO
            ivarp1 = n + 1
  530       CONTINUE
            isofar = ivarp1 - 1
            INVSET( j ) = isofar
          END DO

!  Reorder the elements into set order. Fill in each set from the front. As a
!  new entry is placed in set K increase the pointer ISVSET( k ) by one
!  and find the new variable, INVSET( k ), that corresponds to the set now
!  pointed to by ISVSET( k )

          DO j = 1, nsets + 1

!  Determine the next unplaced entry, ISVSET, in ISET

  560       CONTINUE
            istrt = ISVSET( j )
            IF ( istrt == ISVSET( j + 1 ) ) CYCLE
            IF ( ISET( istrt ) > 0 ) CYCLE

!  Extract the variable and set numbers of the starting element

            ivar = INVSET( j )
            jset = - ISET( istrt )

!  Move elements in a cycle, ending back at set J

            DO k = istrt, n

!  Find the first empty location in set JSET in INVSET

              inext = ISVSET( jset )

!  Extract the variable index of the next element

              newvar = INVSET( jset )

!  Update ISVSET( jset ), find the new variable index and store it in
!  INVSET( jset )

              istarj = inext + 1
              ISVSET( jset ) = istarj
              DO ivarp1 = newvar + 1, n
                 IF ( istarj < ivarp1 ) GO TO 570
              END DO
              ivarp1 = n + 1
  570         CONTINUE
              INVSET( jset ) = ivarp1 - 1
              IF ( jset == j ) EXIT

!  Extract the number of the set of the next element

              newset = - ISET( inext )

!  Store the variable index of the current element

              ISET( inext ) = ivar

!  Make the next element into the current one

              ivar = newvar
              jset = newset
            END DO

!  Store the variable index of the starting element

            ISET( istrt ) = ivar
            GO TO 560
          END DO

!  Revise ISVSET to point to the start of each set

          ISVSET( nsets + 1 : 2 : - 1 ) = ISVSET( nsets : 1 : - 1 )
          ISVSET( 1 ) = 1
        END IF

      ELSE

!  Allocate unused arrays to have length zero

        reallocate = .TRUE.
        IF ( ALLOCATED( ITYPER ) ) THEN
          IF ( SIZE( ITYPER ) /= 0 ) THEN ; DEALLOCATE( ITYPER )
          ELSE ; reallocate = .FALSE. ; END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( ITYPER( 0 ), STAT = alloc_status )
          IF ( alloc_status /= 0 ) THEN
            bad_alloc = 'ITYPER' ; GO TO 600 ; END IF
        END IF

        IF ( ALLOCATED( ISSWTR ) ) THEN
          IF ( SIZE( ISSWTR ) /= 0 ) THEN ; DEALLOCATE( ISSWTR )
          ELSE ; reallocate = .FALSE. ; END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( ISSWTR( 0 ), STAT = alloc_status )
          IF ( alloc_status /= 0 ) THEN
            bad_alloc = 'ISSWTR' ; GO TO 600 ; END IF
        END IF

        reallocate = .TRUE.
        IF ( ALLOCATED( ISSITR ) ) THEN
          IF ( SIZE( ISSITR ) /= 0 ) THEN ; DEALLOCATE( ISSITR )
          ELSE ; reallocate = .FALSE. ; END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( ISSITR( 0 ), STAT = alloc_status )
          IF ( alloc_status /= 0 ) THEN 
            bad_alloc = 'ISSITR' ; GO TO 600 ; END IF
        END IF

        reallocate = .TRUE.
        lwtran = 0
        IF ( ALLOCATED( WTRANS ) ) THEN
          IF ( SIZE( WTRANS ) /= 0 ) THEN ; DEALLOCATE( WTRANS )
          ELSE ; reallocate = .FALSE. ; END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( WTRANS( lwtran ), STAT = alloc_status )
          IF ( alloc_status /= 0 ) THEN
            bad_alloc = 'WTRANS' ; GO TO 600 ; END IF
        END IF
        
        reallocate = .TRUE.
        litran = 0
        IF ( ALLOCATED( ITRANS ) ) THEN
          IF ( SIZE( ITRANS ) /= 0 ) THEN ; DEALLOCATE( ITRANS )
          ELSE ; reallocate = .FALSE. ; END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( ITRANS( litran ), STAT = alloc_status )
          IF ( alloc_status /= 0 ) THEN ; 
            bad_alloc = 'ITRANS' ; GO TO 600 ; END IF
        END IF

        reallocate = .TRUE.
        IF ( ALLOCATED( ISET ) ) THEN
          IF ( SIZE( ISET ) /= 0 ) THEN ; DEALLOCATE( ISET )
          ELSE ; reallocate = .FALSE. ; END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( ISET( 0 ), STAT = alloc_status )
          IF ( alloc_status /= 0 ) THEN 
            bad_alloc = 'ISET' ; GO TO 600 ; END IF
        END IF
        
        reallocate = .TRUE.
        IF ( ALLOCATED( ISVSET ) ) THEN
          IF ( SIZE( ISVSET ) /= 0 ) THEN ; DEALLOCATE( ISVSET )
          ELSE ; reallocate = .FALSE. ; END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( ISVSET( 0 ), STAT = alloc_status )
          IF ( alloc_status /= 0 ) THEN 
            bad_alloc = 'ISVSET' ; GO TO 600 ; END IF
        END IF
        
        reallocate = .TRUE.
        IF ( ALLOCATED( INVSET ) ) THEN
          IF ( SIZE( INVSET ) /= 0 ) THEN ; DEALLOCATE( INVSET )
          ELSE ; reallocate = .FALSE. ; END IF
        END IF
        IF ( reallocate ) THEN 
          ALLOCATE( INVSET( 0 ), STAT = alloc_status )
          IF ( alloc_status /= 0 ) THEN 
            bad_alloc = 'INVSET' ; GO TO 600 ; END IF
        END IF

      END IF

!  Set the length of the remaining partitions of the workspace for array bound
!  checking in calls to other subprograms

!  -- Set the starting addresses for the partitions within FUVALS --

!  A full description of the partitions of FUVALS is given in the introductory
!  comments to the LANCELOT package

      lfxi = 0
      lgxi = lfxi + nel
      lhxi = INTVAR( nel1 ) - 1
      lggfx = lggfx - 1
      ldx = lggfx + n

!  Print all of the starting addresses for the workspace array partitions

       IF ( iprint >= 3 ) WRITE( iout,                                         &
            "( /,' Starting addresses for the partitions of FUVALS ', /,       &
         &       ' ----------------------------------------------- ', //,      &
         &       '   lfxi   lgxi   lhxi  lggfx ','   ldx ', /, 5I7 )" )        &
           lfxi, lgxi, lhxi, lggfx, ldx

!  Set the length of each partition of the real workspace array FUVALS for
!  array bound checking in calls to other subprograms

      lnguvl = MAX( 1, lhxi - lfxi )
      lnhuvl = MAX( 1, lggfx - lfxi )
      status = 0
      RETURN

!  Unsuccessful returns

  600 CONTINUE
      status = 1000 + alloc_status

  610 CONTINUE
      WRITE( iout, 2600 ) bad_alloc, alloc_status
      RETURN

!  Non-executable statements

 2600 FORMAT( ' ** Message from -CUTEST_initialize_workspace-', /,             &
              ' Allocation error, for ', A, ', status = ', I0 )

!  end of subroutine CUTEST_initialize_workspace

      END SUBROUTINE CUTEST_initialize_workspace

!-*-*-*-  C U T E S T _ f o r m _ g r a d i e n t s  S U B R O U T I N E -*-*-*-

      SUBROUTINE CUTEST_form_gradients(                                        &
                       n, ng, nel, ntotel, nvrels, nnza, nvargp,               &
                       firstg, ICNA, ISTADA, IELING, ISTADG, ISTAEV,           &
                       IELVAR, INTVAR, A, GVALS2, GUVALS, lguval,              &
                       GRAD, GSCALE, ESCALE, GRJAC, GXEQX, INTREP,             &
                       ISVGRP, ISTAGV, ITYPEE, ISTAJC, GRAD_el, W_el,          &
                       RANGE, KNDOFG )

!  Calculate the the gradient, GRAD, of the objective function and the
!  Jacobian matrix of gradients, GRJAC, of each group

!  History -
!   ( based on Conn-Gould-Toint fortran 77 version LANCELOT A, ~1992 )
!   fortran 90 version originally released pre GALAHAD Version 1.0. February 
!     7th 1995 as LANCELOT_form_gradients as part of the LANCELOT module
!   update released with GALAHAD Version 2.0. February 16th 2005
!   fortran 2003 version released in CUTEst, 5th November 2012

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER, INTENT( IN    ) :: n, ng, nel, ntotel, nnza, nvargp
      INTEGER, INTENT( IN    ) :: nvrels, lguval
      LOGICAL, INTENT( IN    ) :: firstg
      INTEGER, INTENT( IN    ), DIMENSION( ng  + 1 ) :: ISTADA, ISTADG
      INTEGER, INTENT( IN    ), DIMENSION( nel + 1 ) :: ISTAEV, INTVAR
      INTEGER, INTENT( IN    ), DIMENSION( nvrels  ) :: IELVAR
      INTEGER, INTENT( IN    ), DIMENSION( nnza    ) :: ICNA
      INTEGER, INTENT( IN    ), DIMENSION( ntotel  ) :: IELING
      REAL ( KIND = wp ), INTENT( IN  ), DIMENSION( nnza ) :: A
      REAL ( KIND = wp ), INTENT( IN  ), DIMENSION( ng ) :: GVALS2
      REAL ( KIND = wp ), INTENT( IN  ), DIMENSION( lguval ) :: GUVALS
      REAL ( KIND = wp ), INTENT( IN  ), DIMENSION( ng ) :: GSCALE
      REAL ( KIND = wp ), INTENT( IN  ), DIMENSION( ntotel ) :: ESCALE
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: GRAD
      REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( nvargp ) :: GRJAC
      LOGICAL, INTENT( IN ), DIMENSION( ng  ) :: GXEQX
      LOGICAL, INTENT( IN ), DIMENSION( nel ) :: INTREP
      INTEGER, INTENT( IN ), DIMENSION( : ) :: ISVGRP
      INTEGER, INTENT( IN ), DIMENSION( : ) :: ISTAGV
      INTEGER, INTENT( IN ), DIMENSION( nel ) :: ITYPEE
      INTEGER, INTENT( INOUT ), DIMENSION( : ) :: ISTAJC
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( : ) :: GRAD_el
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( : ) :: W_el
      INTEGER, INTENT( IN ), OPTIONAL, DIMENSION( ng ) :: KNDOFG

!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------

      INTERFACE
        SUBROUTINE RANGE( ielemn, transp, W1, W2, nelvar, ninvar, ieltyp,      &
                          lw1, lw2 )
        INTEGER, INTENT( IN ) :: ielemn, nelvar, ninvar, ieltyp, lw1, lw2
        LOGICAL, INTENT( IN ) :: transp
        REAL ( KIND = KIND( 1.0D+0 ) ), INTENT( IN ), DIMENSION ( lw1 ) :: W1
        REAL ( KIND = KIND( 1.0D+0 ) ), INTENT( OUT ), DIMENSION ( lw2 ) :: W2
        END SUBROUTINE RANGE
      END INTERFACE

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: i, iel, ig, ii, k, ig1, j, jj, l, ll
      INTEGER :: nin, nvarel, nelow, nelup, istrgv, iendgv
      REAL ( KIND = wp ) :: gi, scalee
      LOGICAL :: nontrv

!  Initialize the gradient as zero

      GRAD = zero

!  Consider the IG-th group

      DO ig = 1, ng
        IF ( PRESENT( KNDOFG ) ) THEN
          IF ( KNDOFG( ig ) == 0 ) CYCLE ; END IF
        ig1 = ig + 1
        istrgv = ISTAGV( ig ) ; iendgv = ISTAGV( ig1 ) - 1
        nelow = ISTADG( ig ) ; nelup = ISTADG( ig1 ) - 1
        nontrv = .NOT. GXEQX( ig )

!  Compute the first derivative of the group

        IF ( nontrv ) THEN
          gi = GSCALE( ig ) * GVALS2( ig )
        ELSE
          gi = GSCALE( ig )
        END IF

!  This is the first gradient evaluation or the group has nonlinear elements

        IF ( firstg .OR. nelow <= nelup ) THEN
          GRAD_el( ISVGRP( istrgv : iendgv ) ) = zero

!  Loop over the group's nonlinear elements

          DO ii = nelow, nelup
            iel = IELING( ii )
            k = INTVAR( iel ) ; l = ISTAEV( iel )
            nvarel = ISTAEV( iel + 1 ) - l
            scalee = ESCALE( ii )
            IF ( INTREP( iel ) ) THEN

!  The IEL-th element has an internal representation

              nin = INTVAR( iel + 1 ) - k
              CALL RANGE ( iel, .TRUE., GUVALS( k : k + nin - 1 ),             &
                           W_el( : nvarel ), nvarel, nin, ITYPEE( iel ),       &
                           nin, nvarel )
!DIR$ IVDEP
              DO i = 1, nvarel
                j = IELVAR( l )
                GRAD_el( j ) = GRAD_el( j ) + scalee * W_el( i )
                l = l + 1
              END DO
            ELSE

!  The IEL-th element has no internal representation

!DIR$ IVDEP
              DO i = 1, nvarel
                j = IELVAR( l )
                GRAD_el( j ) = GRAD_el( j ) + scalee * GUVALS( k )
                k = k + 1
                l = l + 1
              END DO
            END IF
          END DO

!  Include the contribution from the linear element

!DIR$ IVDEP
          DO k = ISTADA( ig ), ISTADA( ig1 ) - 1
            GRAD_el( ICNA( k ) ) = GRAD_el( ICNA( k ) ) + A( k )
          END DO

!  Find the gradient of the group

          IF ( nontrv ) THEN

!  The group is non-trivial

!DIR$ IVDEP
            DO i = istrgv, iendgv
              ll = ISVGRP( i )
              GRAD( ll ) = GRAD( ll ) + gi * GRAD_el( ll )

!  As the group is non-trivial, also store the nonzero entries of the
!  gradient of the function in GRJAC

              jj = ISTAJC( ll )
              GRJAC( jj ) = GRAD_el( ll )

!  Increment the address for the next nonzero in the column of
!  the Jacobian for variable LL

              ISTAJC( ll ) = jj + 1
            END DO
          ELSE

!  The group is trivial

!DIR$ IVDEP
            DO i = istrgv, iendgv
              ll = ISVGRP( i )
              GRAD( ll ) = GRAD( ll ) + gi * GRAD_el( ll )
            END DO
          END IF

!  This is not the first gradient evaluation and there is only a linear element

        ELSE

!  Add the gradient of the linear element to the overall gradient

!DIR$ IVDEP
          DO k = ISTADA( ig ), ISTADA( ig1 ) - 1
            GRAD( ICNA( k ) ) = GRAD( ICNA( k ) ) + gi * A( k )
          END DO

!  The group is non-trivial; increment the starting addresses for
!  the groups used by each variable in the (unchanged) linear
!  element to avoid resetting the nonzeros in the Jacobian

          IF ( nontrv ) THEN
!DIR$ IVDEP
            DO i = istrgv, iendgv
              ISTAJC( ISVGRP( i ) ) = ISTAJC( ISVGRP( i ) ) + 1
            END DO
          END IF
        END IF
      END DO

!  Reset the starting addresses for the lists of groups using each variable to
!  their values on entry

      DO i = n, 2, - 1
        ISTAJC( i ) = ISTAJC( i - 1 )
      END DO
      ISTAJC( 1 ) = 1

      RETURN

!  end of subroutine CUTEST_form_gradients

      END SUBROUTINE CUTEST_form_gradients

!-*-*-  C U T E S T _ a s s e m b l e _ h e s s i a n  S U B R O U T I N E -*-*-

      SUBROUTINE CUTEST_assemble_hessian(                                      &
                      n, ng, nel, ntotel, nvrels, nnza, maxsel,                &
                      nvargp, nfree, IFREE, ISTADH, ICNA,                      &
                      ISTADA, INTVAR, IELVAR, IELING, ISTADG, ISTAEV,          &
                      ISTAGV, ISVGRP, A, GUVALS, lnguvl, HUVALS,               &
                      lnhuvl, GVALS2, GVALS3, GSCALE, ESCALE, GXEQX,           &
                      ITYPEE, INTREP, RANGE, iprint, error, out,               &
                      buffer, use_band, no_zeros, fixed_structure, nsemib,     &
                      status, alloc_status, bad_alloc,                         &
                      S, lirnh, ljcnh, lh, IRNH, JCNH, H,                      &
                      LINK_col, POS_in_H, llink, lpos,                         &
                      IW_asmbl, GRAD_el, W_el, W_in, H_el, H_in, skipg,        &
                      nnzh, maxsbw, DIAG, OFFDIA, KNDOFG )

!  Assemble the second derivative matrix of a groups partially separable
!  function in either co-ordinate or band format

!  History -
!   ( based on Conn-Gould-Toint fortran 77 version LANCELOT A, ~1992 )
!   fortran 90 version originally released pre GALAHAD Version 1.0. January
!     25th 1995 as ASMBL_assemble_hessian as part of the ASMBL module
!   update released with GALAHAD Version 2.0. February 16th 2005
!   fortran 2003 version released in CUTEst, 5th November 2012

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

      INTEGER, INTENT( IN  ) :: n, ng, maxsel, nsemib, nvargp, nnza
      INTEGER, INTENT( IN  ) :: nvrels, ntotel, nfree, nel, buffer
      INTEGER, INTENT( IN  ) :: lnguvl, lnhuvl, iprint, error, out
      INTEGER, INTENT( OUT ) :: status, alloc_status
      LOGICAL, INTENT( IN  ) :: use_band, no_zeros, fixed_structure, skipg
      CHARACTER ( LEN = 24 ) :: bad_alloc
      INTEGER, INTENT( IN  ), DIMENSION( n       ) :: IFREE
      INTEGER, INTENT( IN  ), DIMENSION( nnza    ) :: ICNA
      INTEGER, INTENT( IN  ), DIMENSION( ng  + 1 ) :: ISTADA, ISTADG, ISTAGV
      INTEGER, INTENT( IN  ), DIMENSION( nel + 1 ) :: INTVAR, ISTAEV, ISTADH
      INTEGER, INTENT( IN  ), DIMENSION( nvrels  ) :: IELVAR
      INTEGER, INTENT( IN  ), DIMENSION( ntotel  ) :: IELING
      INTEGER, INTENT( IN  ), DIMENSION( nvargp  ) :: ISVGRP
      INTEGER, INTENT( IN  ), DIMENSION( nel ) :: ITYPEE
      REAL ( KIND = wp ), INTENT( IN  ), DIMENSION( nnza ) :: A
      REAL ( KIND = wp ), INTENT( IN  ), DIMENSION( lnguvl ) :: GUVALS
      REAL ( KIND = wp ), INTENT( IN  ), DIMENSION( lnhuvl ) :: HUVALS
      REAL ( KIND = wp ), INTENT( IN  ), DIMENSION( ng ) :: GVALS2
      REAL ( KIND = wp ), INTENT( IN  ), DIMENSION( ng ) :: GVALS3
      REAL ( KIND = wp ), INTENT( IN  ), DIMENSION( ng ) :: GSCALE
      REAL ( KIND = wp ), INTENT( IN  ), DIMENSION( ntotel ) :: ESCALE
      LOGICAL, INTENT( IN ), DIMENSION( ng  ) :: GXEQX
      LOGICAL, INTENT( IN ), DIMENSION( nel ) :: INTREP
      TYPE ( CUTEST_assemble_type ), INTENT( INOUT ) :: S

!---------------------------------------------------------------
!   D u m m y   A r g u m e n t s   f o r   W o r k s p a c e 
!--------------------------------------------------------------

      INTEGER, INTENT( INOUT ) :: lirnh, ljcnh, lh, llink, lpos
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: LINK_col
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: POS_in_H
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: IRNH
      INTEGER, ALLOCATABLE, DIMENSION( : ) :: JCNH
      REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: H

      INTEGER, INTENT( OUT ), DIMENSION( : ) :: IW_asmbl
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( : ) :: GRAD_el
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( : ) :: W_el
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( : ) :: W_in
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( : ) :: H_el
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( : ) :: H_in

!--------------------------------------------------
!   O p t i o n a l   D u m m y   A r g u m e n t s
!--------------------------------------------------

      INTEGER, INTENT( OUT ), OPTIONAL :: maxsbw, nnzh
      REAL ( KIND = wp ), INTENT( OUT ), OPTIONAL,                             &
                                         DIMENSION( nfree ) :: DIAG
      REAL ( KIND = wp ), INTENT( OUT ), OPTIONAL,                             &
                                         DIMENSION( nsemib, nfree ) :: OFFDIA
      INTEGER, INTENT( IN ), OPTIONAL, DIMENSION( ng ) :: KNDOFG

!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------

      INTERFACE
        SUBROUTINE RANGE( ielemn, transp, W1, W2, nelvar, ninvar, ieltyp,      &
                          lw1, lw2 )
        INTEGER, INTENT( IN ) :: ielemn, nelvar, ninvar, ieltyp, lw1, lw2
        LOGICAL, INTENT( IN ) :: transp
        REAL ( KIND = KIND( 1.0D+0 ) ), INTENT( IN  ), DIMENSION ( lw1 ) :: W1
        REAL ( KIND = KIND( 1.0D+0 ) ), INTENT( OUT ), DIMENSION ( lw2 ) :: W2
        END SUBROUTINE RANGE
      END INTERFACE

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: i, ii, j, jj, k, ig, ip, l, newpt, nin
      INTEGER :: iell, iel, ihnext, nvarel, ig1, listvs, listve
      INTEGER :: ielh, inext, ijhess, irow, jcol, jcolst, istart
      INTEGER :: nlh, ulh, mlh
      REAL ( KIND = wp ) :: wki, hesnew, gdash, g2dash, scalee
      CHARACTER ( LEN = 2 ), DIMENSION( 36, 36 ) :: MATRIX
!     CHARACTER ( LEN = 80 ) :: array

!  If a band storage scheme is to be used, initialize the entries within the
!  band as zero

      IF ( use_band ) THEN
        maxsbw = 0
        DIAG = zero ; OFFDIA = zero

!  If a co-ordinate scheme is to be used, allocate arrays to hold the link
!  list which points to the row numbers  which are used in the columns of
!  the assembled Hessian

!  LINK_col( . ) gives the link list. The list for column J starts
!                 in LINK_col( J ) and ends when LINK_col( K ) = - 1
!  POS_in_H( . ) gives the position in H of the current link

      ELSE
        IF ( .NOT. S%ptr_status ) THEN
          S%ptr_status = .TRUE.
          llink = MIN( lirnh, ljcnh, lh ) + nfree
          lpos = MIN( lirnh, ljcnh, lh ) + nfree
        
          ALLOCATE( LINK_col( llink ), STAT = alloc_status )
          IF ( alloc_status /= 0 ) THEN ; bad_alloc = 'LINK_c' ; GO TO 980
          END IF
        
          ALLOCATE( POS_in_H( lpos ), STAT = alloc_status )
          IF ( alloc_status /= 0 ) THEN ; bad_alloc = 'POS_in' ; GO TO 980
          END IF
        
          ALLOCATE( IRNH( lirnh ), STAT = alloc_status )
          IF ( alloc_status /= 0 ) THEN ; bad_alloc = 'IRNH' ; GO TO 980
          END IF
        
          ALLOCATE( JCNH( ljcnh ), STAT = alloc_status )
          IF ( alloc_status /= 0 ) THEN ; bad_alloc = 'JCNH' ; GO TO 980
          END IF
        
          ALLOCATE( H( lh ), STAT = alloc_status )
          IF ( alloc_status /= 0 ) THEN ; bad_alloc = 'H' ; GO TO 980
          END IF
        
        ELSE
        
          k = MIN( lirnh, ljcnh, lh ) + nfree
          IF ( llink < k ) THEN
            DEALLOCATE( LINK_col ) ; llink = k
            ALLOCATE( LINK_col( llink ), STAT = alloc_status )
            IF ( alloc_status /= 0 ) THEN
              bad_alloc = 'LINK_c'; GO TO 980 ; END IF
          END IF
        
          k = MIN( lirnh, ljcnh, lh ) + nfree
          IF ( lpos < k ) THEN 
            DEALLOCATE( POS_in_H ) ; lpos = k
            ALLOCATE( POS_in_H( lpos ), STAT = alloc_status )
            IF ( alloc_status /= 0 ) THEN 
              bad_alloc = 'POS_in'; GO TO 980 ; END IF
          END IF
        
        END IF
        LINK_col( : nfree ) = - 1 ; POS_in_H( : nfree ) = - 1
        newpt = nfree

!  Make an initial allocation of the space required to hold the Hessian

        nnzh = 0
      END IF

!  Renumber the free variables so that they are variables 1 to NFREE

      IW_asmbl( : n ) = 0
      DO i = 1, nfree
        IW_asmbl( IFREE( i ) ) = i
      END DO
      IF ( iprint >= 10 ) WRITE( out, 2060 ) nfree, ( IFREE( i ), i = 1, nfree )

!  ------------------------------------------------------
!  Form the rank-one second order term for the I-th group
!  ------------------------------------------------------

      DO ig = 1, ng
        IF ( skipg ) THEN
          IF ( KNDOFG( ig ) == 0 ) CYCLE ; END IF
        IF ( GXEQX( ig ) ) CYCLE
        IF ( .NOT. fixed_structure .AND. GSCALE( ig ) == zero ) CYCLE
        IF ( iprint >= 100 ) WRITE( out, 2070 ) ig
        g2dash = GSCALE( ig ) * GVALS3( ig )
        IF ( iprint >= 100 ) WRITE( 6, * ) ' GVALS3( ig ) ', GVALS3( ig )
        IF ( no_zeros .AND. g2dash == zero ) CYCLE
        ig1 = ig + 1
        listvs = ISTAGV( ig )
        listve = ISTAGV( ig1 ) - 1

!  Form the gradient of the IG-th group

        GRAD_el( ISVGRP( listvs : listve ) ) = zero

!  Consider any nonlinear elements for the group

        DO iell = ISTADG( ig ), ISTADG( ig1 ) - 1
          iel = IELING( iell )
          k = INTVAR( iel )
          l = ISTAEV( iel )
          nvarel = ISTAEV( iel + 1 ) - l
          scalee = ESCALE( iell )
          IF ( INTREP( iel ) ) THEN

!  The IEL-th element has an internal representation

            nin = INTVAR( iel + 1 ) - k
            CALL RANGE ( iel, .TRUE., GUVALS( k : k + nin - 1 ),               &
                         H_el, nvarel, nin, ITYPEE( iel ), nin, nvarel )
            DO i = 1, nvarel
              j = IELVAR( l )
              GRAD_el( j ) = GRAD_el( j ) + scalee * H_el( i )
              l = l + 1
            END DO
          ELSE

!  The IEL-th element has no internal representation

            DO i = 1, nvarel
              j = IELVAR( l )
              GRAD_el( j ) = GRAD_el( j ) + scalee * GUVALS( k )
              k = k + 1
              l = l + 1
            END DO
          END IF
        END DO

!  Include the contribution from the linear element

        DO k = ISTADA( ig ), ISTADA( ig1 ) - 1
          j = ICNA( k )
          GRAD_el( j ) = GRAD_el( j ) + A( k )
        END DO

!  The gradient is complete. Form the J-TH column of the rank-one matrix

        DO l = listvs, listve
          jj = ISVGRP( l )
          j = IW_asmbl( jj )
          IF ( j == 0 ) CYCLE

!  Find the entry in row I of this column

          DO k = listvs, listve
            ii = ISVGRP( k )
            i = IW_asmbl( ii )
            IF ( i == 0 .OR. i > j ) CYCLE

!  Skip all elements which lie outside a band of width NSEMIB

            IF ( use_band ) maxsbw = MAX( maxsbw, j - i )
            IF ( j - i > nsemib ) CYCLE
            hesnew = GRAD_el( ii ) * GRAD_el( jj ) * g2dash
            IF ( iprint >= 100 ) WRITE( out, 2090 ) i, j, hesnew
            IF ( no_zeros .AND. hesnew == zero ) CYCLE

!  Obtain the appropriate storage location in H for the new entry

!  Case 1: band matrix storage scheme

            IF ( use_band ) THEN

!  The entry belongs on the diagonal

              IF ( i == j ) THEN
                DIAG( i ) = DIAG( i ) + hesnew

!  The entry belongs off the diagonal

              ELSE
                OFFDIA( j - i, i ) = OFFDIA( j - i, i ) + hesnew
              END IF

!  Case 2: co-ordinate storage scheme

            ELSE
              istart = j
 150          CONTINUE
              inext = LINK_col( istart )
              IF ( inext == - 1 ) THEN

!  the (i,j)-th location is empty. Place the new entry in this location and
!  add another link to the list

                nnzh = nnzh + 1
                IF ( nnzh > lh .OR. nnzh > lirnh .OR. nnzh > ljcnh ) THEN
                  nlh = lirnh ; ulh = nnzh - 1 ; mlh = nnzh
                  CALL CUTEST_extend_array( IRNH, lirnh, ulh, nlh, mlh,        &
                                            buffer, status, alloc_status )
                  IF ( status /= 0 ) THEN
                    bad_alloc = 'IRNH' ; GO TO 980 ; END IF
                  lirnh = nlh
                  nlh = ljcnh ; ulh = nnzh - 1 ; mlh = nnzh
                  CALL CUTEST_extend_array( JCNH, lirnh, ulh, nlh, mlh,        &
                                            buffer, status, alloc_status )
                  IF ( status /= 0 ) THEN
                    bad_alloc = 'JCNH' ; GO TO 980 ; END IF
                  ljcnh = nlh
                  nlh = lh ; ulh = nnzh - 1 ; mlh = nnzh
                  CALL CUTEST_extend_array( H, lh, ulh, nlh, mlh, buffer,      &
                                            status, alloc_status )
                  IF ( status /= 0 ) THEN
                    bad_alloc = 'H' ; GO TO 980 ; END IF
                  lh = nlh
                END IF
                IRNH( nnzh ) = i ; JCNH( nnzh ) = j
                H( nnzh ) = hesnew
                IF ( newpt == llink ) THEN
                  nlh = llink
                  ulh = newpt; mlh = llink + 1
                  CALL CUTEST_extend_array( LINK_col, llink, ulh, nlh, mlh,    &
                                            buffer, status, alloc_status )
                  IF ( status /= 0 ) THEN
                    bad_alloc = 'LINK_col' ; GO TO 980 ; END IF
                  llink = nlh
                  nlh = lpos
                  ulh = newpt; mlh = lpos + 1
                  CALL CUTEST_extend_array( POS_in_H, lpos, ulh, nlh, mlh,     &
                                            buffer, status, alloc_status )
                  IF ( status /= 0 ) THEN
                    bad_alloc = 'POS_in_H' ; GO TO 980 ; END IF
                  lpos = nlh
                END IF
                newpt = newpt + 1
                LINK_col( istart ) = newpt
                POS_in_H( istart ) = nnzh
                LINK_col( newpt  ) = - 1
                POS_in_H( newpt ) = - 1
              ELSE

!  Continue searching the linked list for an entry in row I, column J

                IF ( IRNH( POS_in_H( istart ) ) == i ) THEN
                  ip = POS_in_H( istart )
                  H( ip ) = H( ip ) + hesnew
                ELSE
                  istart = inext
                  GO TO 150
                END IF
              END IF
            END IF
          END DO
        END DO
      END DO

!  Reset the workspace array to zero

      W_el( : maxsel ) = zero

! --------------------------------------------------------
! Add on the low rank first order terms for the I-th group
! --------------------------------------------------------

      DO ig = 1, ng
        IF ( skipg ) THEN
          IF ( KNDOFG( ig ) == 0 ) CYCLE ; END IF
        IF ( .NOT. fixed_structure .AND. GSCALE( ig ) == zero ) CYCLE
        IF ( iprint >= 100 ) WRITE( out, 2100 ) ig
        IF ( GXEQX( ig ) ) THEN
          gdash = GSCALE( ig )
        ELSE
          gdash = GSCALE( ig ) * GVALS2( ig )
          IF ( iprint >= 100 ) WRITE( 6, * ) ' GVALS2( ig )', GVALS2( ig )
          IF ( no_zeros .AND. gdash == zero ) CYCLE
        END IF
        ig1 = ig + 1

!  See if the group has any nonlinear elements

        DO iell = ISTADG( ig ), ISTADG( ig1 ) - 1
          iel = IELING( iell )
          listvs = ISTAEV( iel )
          listve = ISTAEV( iel + 1 ) - 1
          nvarel = listve - listvs + 1
          ielh = ISTADH( iel )
          ihnext = ielh
          scalee = ESCALE( iell )
          DO l = listvs, listve
            j = IW_asmbl( IELVAR( l ) )
            IF ( j /= 0 ) THEN

!  The IEL-th element has an internal representation. Compute the J-th column
!  of the element Hessian matrix

              IF ( INTREP( iel ) ) THEN

!  Compute the J-th column of the Hessian

                W_el( l - listvs + 1 ) = one

!  Find the internal variables

                nin = INTVAR( iel + 1 ) - INTVAR( iel )
                CALL RANGE ( iel, .FALSE., W_el, W_in, nvarel, nin,            &
                             ITYPEE( iel ), nvarel, nin )

!  Multiply the internal variables by the element Hessian

                H_in( : nin ) = zero

!  Only the upper triangle of the element Hessian is stored

                  jcolst = ielh - 1
                  DO jcol = 1, nin
                    ijhess = jcolst
                    jcolst = jcolst + jcol
                    wki = W_in( jcol ) * gdash
                    DO irow = 1, nin
                      IF ( irow <= jcol ) THEN
                        ijhess = ijhess + 1
                      ELSE
                        ijhess = ijhess + irow - 1
                      END IF
                      H_in( irow ) = H_in( irow ) + wki * HUVALS( ijhess )
                    END DO
                  END DO

!  Scatter the product back onto the elemental variables

                 CALL RANGE ( iel, .TRUE., H_in, H_el, nvarel, nin,            &
                              ITYPEE( iel ), nin, nvarel )
                 W_el( l - listvs + 1 ) = zero
               END IF

!  Find the entry in row I of this column

               DO k = listvs, l
                 i = IW_asmbl( IELVAR( k ) )

!  Skip all elements which lie outside a band of width nsemib

                 IF ( use_band .AND. i /= 0 )                                  &
                   maxsbw = MAX( maxsbw, ABS( j - i ) )
                 IF ( ABS( i - j ) <= nsemib .AND. i /= 0 ) THEN

!  Only the upper triangle of the matrix is stored

                   IF ( i <= j ) THEN
                     ii = i
                     jj = j
                   ELSE
                     ii = j
                     jj = i
                   END IF

!  Obtain the appropriate storage location in H for the new entry

                   IF ( INTREP( iel ) ) THEN
                     hesnew = scalee * H_el( k - listvs + 1 )
                   ELSE
                     hesnew = scalee * HUVALS( ihnext ) * gdash
                   END IF
                   IF ( iprint >= 100 ) WRITE( 6, 2080 ) ii, jj, iel, hesnew

!  Case 1: band matrix storage scheme

                   IF ( use_band ) THEN

!  The entry belongs on the diagonal

                     IF ( ii == jj ) THEN
                       DIAG( ii ) = DIAG( ii ) + hesnew
                       IF ( k /= l ) DIAG( ii ) = DIAG( ii ) + hesnew

!  The entry belongs off the diagonal

                     ELSE
                       OFFDIA( jj - ii, ii ) = OFFDIA( jj - ii, ii ) + hesnew
                     END IF

!  Case 2: co-ordinate storage scheme

                   ELSE
                     IF ( .NOT. no_zeros .OR. hesnew /= zero ) THEN
                      istart = jj
  230                 CONTINUE
                      inext = LINK_col( istart )
                      IF ( inext == - 1 ) THEN

!  The (i,j)-th location is empty. Place the new entry in this location
!  and add another link to the list

                        nnzh = nnzh + 1
                        
                        IF ( nnzh > lh .OR.                                    &
                             nnzh > lirnh .OR. nnzh > ljcnh ) THEN
                          nlh = lirnh ; ulh = nnzh - 1; mlh = nnzh
                          CALL CUTEST_extend_array( IRNH, lirnh, ulh, nlh,     &
                                                    mlh, buffer, status,       &
                                                    alloc_status )
                          IF ( status /= 0 ) THEN
                             bad_alloc = 'IRNH' ; GO TO 980 ; END IF
                          lirnh = nlh
                          nlh = ljcnh ; ulh = nnzh - 1; mlh = nnzh
                          CALL CUTEST_extend_array( JCNH, lirnh, ulh, nlh,     &
                                                    mlh, buffer, status,       &
                                                    alloc_status )
                          IF ( status /= 0 ) THEN
                             bad_alloc = 'JCNH' ; GO TO 980 ; END IF
                          ljcnh = nlh
                          nlh = lh ; ulh = nnzh - 1 ; mlh = nnzh
                          CALL CUTEST_extend_array( H, lh, ulh, nlh,           &
                                                    mlh, buffer, status,       &
                                                    alloc_status )
                          IF ( status /= 0 ) THEN
                             bad_alloc = 'H' ; GO TO 980 ; END IF
                          lh = nlh
                        END IF   
                        IRNH( nnzh ) = ii ; JCNH( nnzh ) = jj
                        H( nnzh ) = hesnew
                        IF( k /= l .AND. ii == jj ) H( nnzh ) = hesnew + hesnew
                        IF ( newpt == llink ) THEN
                          nlh = llink
                          ulh = newpt; mlh = llink + 1
                          CALL CUTEST_extend_array( LINK_col, llink, ulh, nlh, &
                                                    mlh, buffer, status,       &
                                                    alloc_status )
                          IF ( status /= 0 ) THEN
                            bad_alloc = 'LINK_col' ; GO TO 980 ; END IF
                          llink = nlh
                          nlh = lpos
                          ulh = newpt; mlh = lpos + 1
                          CALL CUTEST_extend_array( POS_in_H, lpos, ulh, nlh,  &
                                                    mlh, buffer, status,       &
                                                    alloc_status )
                          IF ( status /= 0 ) THEN
                            bad_alloc = 'POS_in_H' ; GO TO 980 ; END IF
                          lpos = nlh
                        END IF
                        newpt = newpt + 1
                        LINK_col( istart ) = newpt
                        POS_in_H( istart ) = nnzh
                        LINK_col( newpt  ) = - 1
                        POS_in_H( newpt ) = - 1
                      ELSE

! Continue searching the linked list for an entry in row I, column J

                        IF ( IRNH( POS_in_H( istart ) ) == ii ) THEN
                          ip = POS_in_H( istart )
                          H( ip ) = H( ip ) + hesnew
                          IF( k /= l .AND. ii == jj ) H( ip ) = H( ip ) + hesnew
                        ELSE
                          istart = inext
                          GO TO 230
                        END IF
                      END IF
                    END IF
                  END IF
                END IF
                ihnext = ihnext + 1
              END DO
            ELSE
              ihnext = ihnext + l - listvs + 1
            END IF
          END DO
        END DO
      END DO

!  ---------------------------------------
!  For debugging, print the nonzero values
!  ---------------------------------------

      IF ( iprint >= 10 ) THEN
        IF ( .NOT. use_band )                                                  &
          WRITE( out, 2000 ) ( IRNH( i ), JCNH( i ), H( i ), i = 1, nnzh )

!  For debugging, print the nonzero pattern of the matrix

        IF ( nfree <= 36 ) THEN
          MATRIX( : nfree, : nfree ) = '  '
          IF ( use_band ) THEN
            DO i = 1, nfree
              IF ( DIAG( i ) /= zero ) MATRIX( i, i ) = ' *'
              DO j = 1, MIN( nsemib, nfree - i )
                IF ( OFFDIA( j, i ) /= zero ) THEN
                   MATRIX( i + j, i ) = ' *'
                   MATRIX( i, i + j ) = ' *'
                END IF
              END DO
            END DO
          ELSE
            DO i = 1, nnzh
              IF ( IRNH( i ) > nfree ) THEN
                WRITE( out, * ) ' Entry out of bounds in ASMBL ',              &
                                 ' row number = ', IRNH( i )
!               STOP
              END IF
              IF ( JCNH( i ) > nfree ) THEN
                WRITE( out, * ) ' Entry out of bounds in ASMBL ',              &
                                 ' col number = ', JCNH( i )
!               STOP
              END IF
              MATRIX( IRNH( i ), JCNH( i ) ) = ' *'
              MATRIX( JCNH( i ), IRNH( i ) ) = ' *'
            END DO
          END IF
          WRITE( out, 2040 ) ( i, i = 1, nfree )
          DO i = 1, nfree
            WRITE( out, 2050 ) i, ( MATRIX( i, j ), j = 1, nfree )
          END DO
        END IF
      END IF

!  Successful return

      status = 0
      RETURN

!  Unsuccessful returns

  980 CONTINUE
      WRITE( error, 2990 ) alloc_status, bad_alloc
      RETURN

!  Non-executable statements

 2000 FORMAT( '    Row  Column    Value        Row  Column    Value ', /       &
              '    ---  ------    -----        ---  ------    ----- ', /       &
              ( 2I6, ES24.16, 2I6, ES24.16 ) )
 2040 FORMAT( /, 5X, 36I2 )
 2050 FORMAT( I3, 2X, 36A2 )
 2060 FORMAT( /, I6, ' free variables. They are ', 8I5, /, ( 14I5 ) )
 2070 FORMAT( ' Group ', I5, ' rank-one terms ' )
 2080 FORMAT( ' Row ', I6, ' Column ', I6, ' used from element ', I6,          &
              ' value = ', ES24.16 )
 2090 FORMAT( ' Row ', I6, ' column ', I6, ' used. Value = ', ES24.16 )
 2100 FORMAT( ' Group ', I5, ' second-order terms ' )
 2990 FORMAT( ' ** Message from -CUTEST_assemble_hessian-', /,                 &
              ' Allocation error (status = ', I6, ') for ', A24 )

!  end of subroutine CUTEST_assemble_hessian

     END SUBROUTINE CUTEST_assemble_hessian

!-  C U T E S T _ h e s s i a n _ t i m e s _ v e c t o r  S U B R O U T I N E -

     SUBROUTINE CUTEST_hessian_times_vector(                                   &
                      n, ng, nel, ntotel, nvrels, nvargp,  nfree, nvar1,       &
                      nvar2, nnonnz, nbprod, alllin, IVAR, ISTAEV, ISTADH,     &
                      INTVAR, IELING, IELVAR, ISWKSP, INONNZ, P, Q, GVALS2,    &
                      GVALS3, GRJAC, GSCALE, ESCALE, HUVALS, lhuval, GXEQX,    &
                      INTREP, densep, IGCOLJ, ISLGRP, ISVGRP, ISTAGV, IVALJR,  &
                      ITYPEE, ISYMMH, ISTAJC, IUSED, LIST_elements,            &
                      LINK_elem_uses_var, NZ_comp_w, AP, W_el, W_in, H_in,     &
                      RANGE, skipg, KNDOFG )

!  ---------------------------------------------------------------------------
!  evaluate Q, the product of the hessian of a groups partially separable
!  function with the vector P

!  the nonzero components of P have indices IVAR(i), i = nvar1, ..., nvar2.
!  The nonzero components of the product Q have indices INNONZ(i),
!  i = 1, ..., nnonnz

!  the elements of the array IUSED must be set to zero on entry; they will have
!  been reset to zero on exit. 

!  the components of ISWKSP must be less than nbprod on entry; on exit they 
!  will be no larger than nbprod
!  ---------------------------------------------------------------------------

!  History -
!   fortran 77 version originally released in CUTE, September 23rd, 1991
!   fortran 90 version originally released pre GALAHAD Version 1.0. Febrauary
!     1st 1995 as HSPRD_hessian_times_vector as part of the HSPRD module
!   update released with GALAHAD Version 2.0. February 16th 2005
!   fortran 2003 version released in CUTEst, 5th November 2012

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

     INTEGER, INTENT( IN    ) :: n, ng, nel, ntotel, nvrels, nfree
     INTEGER, INTENT( IN    ) :: nvar1, nvar2, nbprod
     INTEGER, INTENT( IN    ) :: nvargp, lhuval
     INTEGER, INTENT( INOUT ) :: nnonnz
     LOGICAL, INTENT( IN    ) :: alllin, densep, skipg
     INTEGER, INTENT( IN    ), DIMENSION( n ) :: IVAR
     INTEGER, INTENT( IN    ), DIMENSION( nel + 1 ) :: ISTAEV, ISTADH
     INTEGER, INTENT( IN    ), DIMENSION( nel + 1 ) :: INTVAR
     INTEGER, INTENT( IN    ), DIMENSION( ntotel  ) :: IELING
     INTEGER, INTENT( IN    ), DIMENSION( nvrels  ) :: IELVAR
     INTEGER, INTENT( IN    ), DIMENSION( nel     ) :: ITYPEE
     INTEGER, INTENT( INOUT ), DIMENSION( ntotel ) :: ISWKSP
     INTEGER, INTENT( INOUT ), DIMENSION( n ) :: INONNZ
     REAL ( KIND = wp ), INTENT( IN  ), DIMENSION( n ) :: P
     REAL ( KIND = wp ), INTENT( IN  ), DIMENSION( ng ) :: GVALS2
     REAL ( KIND = wp ), INTENT( IN  ), DIMENSION( ng ) :: GVALS3
     REAL ( KIND = wp ), INTENT( IN  ), DIMENSION( ng ) :: GSCALE
     REAL ( KIND = wp ), INTENT( IN  ), DIMENSION( nvargp ) :: GRJAC
     REAL ( KIND = wp ), INTENT( IN  ), DIMENSION( ntotel ) :: ESCALE
     REAL ( KIND = wp ), INTENT( IN  ), DIMENSION( lhuval ) :: HUVALS
     REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: Q
     LOGICAL, INTENT( IN ), DIMENSION( ng ) :: GXEQX
     LOGICAL, INTENT( IN ), DIMENSION( nel ) :: INTREP
     INTEGER, INTENT( IN ), DIMENSION( : ) :: IGCOLJ
     INTEGER, INTENT( IN ), DIMENSION( : ) :: ISLGRP
     INTEGER, INTENT( IN ), DIMENSION( : ) :: ISVGRP
     INTEGER, INTENT( IN ), DIMENSION( : ) :: ISTAGV
     INTEGER, INTENT( IN ), DIMENSION( : ) :: IVALJR
     INTEGER, INTENT( IN ), DIMENSION( : ) :: ISTAJC
     INTEGER, INTENT( INOUT ), DIMENSION( : ) :: IUSED 
     INTEGER, INTENT( IN ), DIMENSION( : ) :: LIST_elements
     INTEGER, INTENT( IN ), DIMENSION( : , : ) :: ISYMMH

     INTEGER, INTENT( IN ), DIMENSION( : ) :: LINK_elem_uses_var
     INTEGER, INTENT( OUT ), DIMENSION( : ) :: NZ_comp_w
     REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( : ) :: AP
     REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( : ) :: W_el
     REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( : ) :: W_in
     REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( : ) :: H_in

     INTEGER, INTENT( IN ), OPTIONAL, DIMENSION( ng ) :: KNDOFG

!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------

     INTERFACE
       SUBROUTINE RANGE( ielemn, transp, W1, W2, nelvar, ninvar, ieltyp,       &
                         lw1, lw2 )
       INTEGER, INTENT( IN ) :: ielemn, nelvar, ninvar, ieltyp, lw1, lw2
       LOGICAL, INTENT( IN ) :: transp
       REAL ( KIND = KIND( 1.0D+0 ) ), INTENT( IN  ), DIMENSION ( lw1 ) :: W1
       REAL ( KIND = KIND( 1.0D+0 ) ), INTENT( OUT ), DIMENSION ( lw2 ) :: W2
       END SUBROUTINE RANGE
     END INTERFACE

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

     INTEGER :: i, iel, ig, ii, ipt, j, irow, jcol, ijhess, lthvar
     INTEGER :: iell, nin, k, l, ll, nvarel, ielhst, nnz_comp_w
     REAL ( KIND = wp ) :: pi, gi, smallest
     LOGICAL :: nullwk

     smallest = TINY( one )

!  ======================= rank-one terms ==========================

!  If the IG-th group is non-trivial, form the product of P with the
!  sum of rank-one first order terms, A(trans) * GVALS3 * A. A is
!  stored by both rows and columns. For maximum efficiency, the
!  product is formed in different ways if P is sparse or dense

!  -----------------  Case 1. P is not sparse -----------------------

     IF ( densep ) THEN

!  Initialize AP and Q as zero

       AP( : ng ) = zero ; Q = zero

!  Form the matrix-vector product AP = A * P, using the column-wise
!  storage of A

       DO j = nvar1, nvar2
         i = IVAR( j )
         pi = P( i )
!DIR$ IVDEP
         DO k = ISTAJC( i ), ISTAJC( i + 1 ) - 1
           AP( IGCOLJ( k ) ) = AP( IGCOLJ( k ) ) + pi * GRJAC( k )
         END DO
       END DO

!  Multiply W by the diagonal matrix GVALS3

       IF ( skipg ) THEN
         DO ig = 1, ng
           IF ( KNDOFG( ig ) == 0 ) THEN
             AP( ig ) = zero
           ELSE
             IF ( GXEQX( ig ) ) THEN
               AP( ig ) = AP( ig ) * GSCALE( ig )
             ELSE
               AP( ig ) = AP( ig ) * GSCALE( ig ) * GVALS3( ig )
             END IF
           END IF
         END DO
       ELSE
         WHERE ( GXEQX( : ng ) ) ; AP( : ng ) = AP( : ng ) * GSCALE( : ng )
         ELSEWHERE ; AP( : ng ) = AP( : ng ) * GSCALE( : ng ) * GVALS3( : ng )
         END WHERE
       END IF

!  Form the matrix-vector product Q = A(trans) * W, once again using the
!  column-wise storage of A

       nnonnz = 0
       DO j = 1, nfree
         i = IVAR( j )
!        Q( i ) =                                                              &
!          DOT_PRODUCT( AP( IGCOLJ( ISTAJC( i ) : ISTAJC( i + 1 ) - 1 ) ),     &
!                       GRJAC ( ISTAJC( i ) : ISTAJC( i + 1 ) - 1 ) )
         pi = zero
         DO ii = ISTAJC( i ), ISTAJC( i + 1 ) - 1
           pi = pi + AP( IGCOLJ( ii ) ) * GRJAC( ii )
         END DO
         Q( i ) = pi
       END DO
     ELSE

!  ------------------- Case 2. P is sparse --------------------------

       nnz_comp_w = 0
       Q( IVAR( : nfree ) ) = zero

!  Form the matrix-vector product W = A * P, using the column-wise
!  storage of A. Keep track of the nonzero components of W in NZ_comp_w.
!  Only store components corresponding to non trivial groups

       DO j = nvar1, nvar2
         i = IVAR( j )
         pi = P( i )
!DIR$ IVDEP
         DO k = ISTAJC( i ), ISTAJC( i + 1 ) - 1
           ig = IGCOLJ( k )
           IF ( IUSED( ig ) == 0 ) THEN
             AP( ig ) = pi * GRJAC( k )
             IUSED( ig ) = 1
             nnz_comp_w = nnz_comp_w + 1
             NZ_comp_w( nnz_comp_w ) = ig
           ELSE
             AP( ig ) = AP( ig ) + pi * GRJAC( k )
           END IF
         END DO
       END DO

!  Reset IUSED to zero

       IUSED( NZ_comp_w( : nnz_comp_w ) ) = 0

!  Form the matrix-vector product Q = A( TRANS ) * W, using the row-wise
!  storage of A

       nnonnz = 0
       IF ( skipg ) THEN
         DO j = 1, nnz_comp_w
           ig = NZ_comp_w( j )
           IF ( KNDOFG( ig ) == 0 ) CYCLE
           IF ( .NOT. GXEQX( ig ) ) THEN

!  If group ig is non trivial, there are contributions from its rank-one term

             pi = GSCALE( ig ) * GVALS3( ig ) * AP( ig )
!DIR$ IVDEP
             DO k = ISTAGV( ig ), ISTAGV( ig + 1 ) - 1
               l = ISVGRP( k )

!  If Q has a nonzero in position L, store its index in INONNZ

               IF ( IUSED( l ) == 0 ) THEN
                 Q( l ) = pi * GRJAC( IVALJR( k ) )
                 IUSED( l ) = 1
                 nnonnz = nnonnz + 1
                 INONNZ( nnonnz ) = l
               ELSE
                 Q( l ) = Q( l ) + pi * GRJAC( IVALJR( k ) )
               END IF
             END DO
           END IF
         END DO
       ELSE
         DO j = 1, nnz_comp_w
           ig = NZ_comp_w( j )
           IF ( .NOT. GXEQX( ig ) ) THEN

!  If group ig is non trivial, there are contributions from its rank-one term

             pi = GSCALE( ig ) * GVALS3( ig ) * AP( ig )
!DIR$ IVDEP
             DO k = ISTAGV( ig ), ISTAGV( ig + 1 ) - 1
               l = ISVGRP( k )

!  If Q has a nonzero in position L, store its index in INONNZ

               IF ( IUSED( l ) == 0 ) THEN
                 Q( l ) = pi * GRJAC( IVALJR( k ) )
                 IUSED( l ) = 1
                 nnonnz = nnonnz + 1
                 INONNZ( nnonnz ) = l
               ELSE
                 Q( l ) = Q( l ) + pi * GRJAC( IVALJR( k ) )
               END IF
             END DO
           END IF
         END DO
       END IF
     END IF

     IF ( .NOT. alllin ) THEN

!  ======================= second-order terms =======================

!  Now consider the product of P with the second order terms (that is, the
!  2nd derivatives of the elements). Again, for maximum efficiency, the
!  product is formed in different ways if P is sparse or dense

!  --------------------- Case 1. P is not sparse ---------------------

       IF ( densep ) THEN
         DO iell = 1, ntotel
           ig = ISLGRP( iell )
           IF ( skipg ) THEN ; IF ( KNDOFG( ig ) == 0 ) CYCLE ; END IF
           ISWKSP( iell ) = nbprod
           iel = IELING( iell )
           nvarel = ISTAEV( iel + 1 ) - ISTAEV( iel )
           IF ( GXEQX( ig ) ) THEN
             gi = GSCALE( ig ) * ESCALE( iell )
           ELSE
             gi = GSCALE( ig ) * ESCALE( iell ) * GVALS2( ig )
           END IF
           IF ( INTREP( iel ) ) THEN

!  The IEL-th element Hessian has an internal representation. Copy the
!  elemental variables into W

             nullwk = .TRUE.
             ll = ISTAEV( iel )
!DIR$ IVDEP
             DO ii = 1, nvarel
               pi = P( IELVAR( ll ) )
               W_el( ii ) = pi
               IF ( pi /= zero ) nullwk = .FALSE.
               ll = ll + 1
             END DO
             IF ( nullwk ) CYCLE

!  Find the internal variables, W_in

             nin = INTVAR( iel + 1 ) - INTVAR( iel )
             CALL RANGE ( iel, .FALSE., W_el, W_in, nvarel, nin,               &
                          ITYPEE( iel ), nvarel, nin )

!  Multiply the internal variables by the element Hessian and put the
!  product in H_in. Consider the first column of the element Hessian

             ielhst = ISTADH( iel )
             pi = gi * W_in( 1 )
             H_in( : nin ) = pi * HUVALS( ISYMMH( 1, : nin ) + ielhst )

!  Now consider the remaining columns of the element Hessian

             DO jcol = 2, nin
               pi = gi * W_in( jcol )
               IF ( pi /= zero ) THEN
                 H_in( : nin ) = H_in( : nin ) +                               &
                   pi * HUVALS( ISYMMH( jcol, : nin )+ ielhst )
               END IF
             END DO

!  Scatter the product back onto the elemental variables, W

             CALL RANGE ( iel, .TRUE., H_in, W_el, nvarel, nin,                &
                          ITYPEE( iel ), nin, nvarel )

!  Add the scattered product to Q

             ll = ISTAEV( iel )
!DIR$ IVDEP
             DO ii = 1, nvarel
                l = IELVAR( ll )
                Q( l ) = Q( l ) + W_el( ii )
                ll = ll + 1
             END DO
           ELSE

!  The IEL-th element Hessian has no internal representation

             lthvar = ISTAEV( iel ) - 1
             ielhst = ISTADH( iel )
             DO jcol = 1, nvarel
               pi = gi * P( IELVAR( lthvar + jcol ) )
               IF ( pi /= zero ) THEN
!DIR$ IVDEP  
                 DO irow = 1, nvarel
                   ijhess = ISYMMH( jcol, irow ) + ielhst
                   l = IELVAR( lthvar + irow )
                   Q( l ) = Q( l ) + pi * HUVALS( ijhess )
                 END DO
               END IF
             END DO
           END IF
         END DO
       ELSE

!  -------------------- Case 2. P is sparse ------------------------

         IF ( skipg ) THEN 
           DO j = nvar1, nvar2

!  Consider each nonzero component of P separately

             i = IVAR( j )
             ipt = LINK_elem_uses_var( i )
             IF ( ipt >= 0 ) THEN

!  The index of the I-th component lies in the IEL-th nonlinear element

               iell = LIST_elements( i )
  300          CONTINUE

!  Check to ensure that the IEL-th element has not already been used

               IF ( ISWKSP( iell ) < nbprod ) THEN
                 ISWKSP( iell ) = nbprod
                 ig = ISLGRP( iell )
                 IF ( KNDOFG( ig ) /= 0 ) THEN
                   iel = IELING( iell )
                   nvarel = ISTAEV( iel + 1 ) - ISTAEV( iel )
                   IF ( GXEQX( ig ) ) THEN
                     gi = GSCALE( ig ) * ESCALE( iell )
                   ELSE
                     gi = GSCALE( ig ) * ESCALE( iell ) * GVALS2( ig )
                   END IF
                   IF ( INTREP( iel ) ) THEN

!  The IEL-th element Hessian has an internal representation. Copy the
!  elemental variables into W

                     ll = ISTAEV( iel )
                     W_el( : nvarel ) = P( IELVAR( ll : ll + nvarel - 1 ) )

!  Find the internal variables

                     nin = INTVAR( iel + 1 ) - INTVAR( iel )
                     CALL RANGE ( iel, .FALSE., W_el, W_in, nvarel, nin,       &
                                  ITYPEE( iel ), nvarel, nin )

!  Multiply the internal variables by the element Hessian and put the
!  product in W_in. Consider the first column of the element Hessian

                     ielhst = ISTADH( iel )
                     pi = gi * W_in( 1 )
                     H_in( : nin ) = pi * HUVALS( ISYMMH( 1, : nin ) + ielhst )

!  Now consider the remaining columns of the element Hessian

                     DO jcol = 2, nin
                       pi = gi * W_in( jcol )
                       IF ( pi /= zero ) THEN
                         H_in( : nin ) = H_in( : nin ) + pi *                  &
                           HUVALS( ISYMMH( jcol, : nin ) + ielhst )
                       END IF
                     END DO

!  Scatter the product back onto the elemental variables, W

                     CALL RANGE ( iel, .TRUE., H_in, W_el, nvarel, nin,        &
                                  ITYPEE( iel ), nin, nvarel )

!  Add the scattered product to Q

                     ll = ISTAEV( iel )
!DIR$ IVDEP
                     DO ii = 1, nvarel
                       l = IELVAR( ll )

!  If Q has a nonzero in position L, store its index in INONNZ

                       IF ( ABS( W_el( ii ) ) > smallest ) THEN
                         IF ( IUSED( l ) == 0 ) THEN
                           Q( l ) = W_el( ii )
                           IUSED( l ) = 1
                           nnonnz = nnonnz + 1
                           INONNZ( nnonnz ) = l
                         ELSE
                           Q( l ) = Q( l ) + W_el( ii )
                         END IF
                       END IF
                       ll = ll + 1
                     END DO
                   ELSE

!  The IEL-th element Hessian has no internal representation

                     lthvar = ISTAEV( iel ) - 1
                     ielhst = ISTADH( iel )
                     DO jcol = 1, nvarel
                       pi = gi * P( IELVAR( lthvar + jcol ) )
                       IF ( pi /= zero ) THEN
!DIR$ IVDEP        
                         DO irow = 1, nvarel
                           ijhess = ISYMMH( jcol, irow ) + ielhst

!  If Q has a nonzero in position L, store its index in INONNZ

                           IF ( ABS( HUVALS( ijhess ) ) > smallest ) THEN
                             l = IELVAR( lthvar + irow )
                             IF ( IUSED( l ) == 0 ) THEN
                               Q( l ) = pi * HUVALS( ijhess )
                               IUSED( l ) = 1
                               nnonnz = nnonnz + 1
                               INONNZ( nnonnz ) = l
                             ELSE
                                Q( l ) = Q( l ) + pi * HUVALS( ijhess )
                             END IF
                           END IF
                         END DO
                       END IF
                     END DO
                   END IF
                 END IF
               END IF

!  Check to see if there are any further elements whose variables
!  include the I-th variable

               IF ( ipt > 0 ) THEN
                 iell = LIST_elements( ipt )
                 ipt = LINK_elem_uses_var( ipt )
                 GO TO 300
               END IF
             END IF
           END DO
         ELSE
           DO j = nvar1, nvar2

!  Consider each nonzero component of P separately

             i = IVAR( j )
             ipt = LINK_elem_uses_var( i )
             IF ( ipt >= 0 ) THEN

!  The index of the I-th component lies in the IEL-th nonlinear element

               iell = LIST_elements( i )
  310          CONTINUE

!  Check to ensure that the IEL-th element has not already been used

               IF ( ISWKSP( iell ) < nbprod ) THEN
                 ISWKSP( iell ) = nbprod
                 iel = IELING( iell )
                 nvarel = ISTAEV( iel + 1 ) - ISTAEV( iel )
                 ig = ISLGRP( iell )
                 IF ( GXEQX( ig ) ) THEN
                   gi = GSCALE( ig ) * ESCALE( iell )
                 ELSE
                   gi = GSCALE( ig ) * ESCALE( iell ) * GVALS2( ig )
                 END IF
                 IF ( INTREP( iel ) ) THEN

!  The IEL-th element Hessian has an internal representation. Copy the
!  elemental variables into W

                   ll = ISTAEV( iel )
                   W_el( : nvarel ) = P( IELVAR( ll : ll + nvarel - 1 ) )

!  Find the internal variables

                   nin = INTVAR( iel + 1 ) - INTVAR( iel )
                   CALL RANGE ( iel, .FALSE., W_el, W_in, nvarel, nin,         &
                                ITYPEE( iel ), nvarel, nin )

!  Multiply the internal variables by the element Hessian and put the
!  product in W_in. Consider the first column of the element Hessian

                   ielhst = ISTADH( iel )
                   pi = gi * W_in( 1 )
                   H_in( : nin ) = pi * HUVALS( ISYMMH( 1, : nin ) + ielhst )

!  Now consider the remaining columns of the element Hessian

                   DO jcol = 2, nin
                     pi = gi * W_in( jcol )
                     IF ( pi /= zero ) THEN
                       H_in( : nin ) = H_in( : nin ) + pi *                    &
                         HUVALS( ISYMMH( jcol, : nin ) + ielhst )
                     END IF
                   END DO

!  Scatter the product back onto the elemental variables, W

                   CALL RANGE ( iel, .TRUE., H_in, W_el, nvarel, nin,          &
                                ITYPEE( iel ), nin, nvarel )

!  Add the scattered product to Q

                   ll = ISTAEV( iel )
!DIR$ IVDEP
                   DO ii = 1, nvarel
                     l = IELVAR( ll )

!  If Q has a nonzero in position L, store its index in INONNZ

                     IF ( ABS( W_el( ii ) ) > smallest ) THEN
                       IF ( IUSED( l ) == 0 ) THEN
                         Q( l ) = W_el( ii )
                         IUSED( l ) = 1
                         nnonnz = nnonnz + 1
                         INONNZ( nnonnz ) = l
                       ELSE
                         Q( l ) = Q( l ) + W_el( ii )
                       END IF
                     END IF
                     ll = ll + 1
                   END DO
                 ELSE

!  The IEL-th element Hessian has no internal representation

                   lthvar = ISTAEV( iel ) - 1
                   ielhst = ISTADH( iel )
                   DO jcol = 1, nvarel
                     pi = gi * P( IELVAR( lthvar + jcol ) )
                     IF ( pi /= zero ) THEN
!DIR$ IVDEP      
                       DO irow = 1, nvarel
                         ijhess = ISYMMH( jcol, irow ) + ielhst

!  If Q has a nonzero in position L, store its index in INONNZ

                         IF ( ABS( HUVALS( ijhess ) ) > smallest ) THEN
                           l = IELVAR( lthvar + irow )
                           IF ( IUSED( l ) == 0 ) THEN
                             Q( l ) = pi * HUVALS( ijhess )
                             IUSED( l ) = 1
                             nnonnz = nnonnz + 1
                             INONNZ( nnonnz ) = l
                           ELSE
                              Q( l ) = Q( l ) + pi * HUVALS( ijhess )
                           END IF
                         END IF
                       END DO
                     END IF
                   END DO
                 END IF
               END IF

!  Check to see if there are any further elements whose variables
!  include the I-th variable

               IF ( ipt > 0 ) THEN
                 iell = LIST_elements( ipt )
                 ipt = LINK_elem_uses_var( ipt )
                 GO TO 310
               END IF
             END IF
           END DO
         END IF
       END IF
     END IF

!  ==================== the product is complete =======================

!  Reset IUSED to zero

     IF ( .NOT. densep ) IUSED( INONNZ( : nnonnz ) ) = 0
     RETURN

!  end of subroutine CUTEST_hessian_times_vector

     END SUBROUTINE CUTEST_hessian_times_vector

!-*-  C U T E S T _ e x t e n d _ a r r a y _ r e a l  S U B R O U T I N E -*-

     SUBROUTINE CUTEST_extend_array_real( ARRAY, old_length, used_length,      &
                                          new_length, min_length, buffer,      &
                                          status, alloc_status )

!  ---------------------------------------------------------------------
!  extend a real array so that its length is increaed from old_length to 
!  as close to new_length as possible while keeping existing data intact
!  ---------------------------------------------------------------------

!  History -
!   fortran 90 version released pre GALAHAD Version 1.0. February 7th 1995 as
!     EXTEND_array_real as part of the GALAHAD module EXTEND
!   fortran 2003 version released in CUTEst, 5th November 2012

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

     INTEGER, INTENT( IN ) :: old_length, buffer
     INTEGER, INTENT( OUT ) :: status, alloc_status
     INTEGER, INTENT( INOUT ) :: used_length, min_length, new_length
     REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: ARRAY

     INTEGER :: length
     LOGICAL :: file_open
     REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DUMMY

!  Make sure that the new length is larger than the old

     IF ( new_length <= old_length ) new_length = 2 * old_length

!  Ensure that the input data is consistent

     used_length = MIN( used_length, old_length )
     min_length = MAX( old_length + 1, MIN( min_length, new_length ) )

!  If possible, allocate DUMMY to hold the old values of ARRAY

     ALLOCATE( DUMMY( used_length ), STAT = alloc_status )

!  If the allocation failed, resort to using an external unit

     IF ( alloc_status /= 0 ) GO TO 100

     DUMMY( : used_length ) = ARRAY( : used_length )

!  Extend the length of ARRAY

     DEALLOCATE( ARRAY )
     length = new_length

  10 CONTINUE
     ALLOCATE( ARRAY( length ), STAT = alloc_status )

!  If the allocation failed, reduce the new length and retry

     IF ( alloc_status /= 0 ) THEN
       length = length + ( length - min_length ) / 2

!  If there is insufficient room for both ARRAY and DUMMY, use an external unit

       IF ( length < min_length ) THEN

!  Rewind the buffer i/o unit

         INQUIRE( UNIT = buffer, OPENED = file_open )
         IF ( file_open ) THEN
           REWIND( UNIT = buffer )
         ELSE
           OPEN( UNIT = buffer )
         END IF

!  Copy the contents of ARRAY into the buffer i/o area

         WRITE( UNIT = buffer, FMT = * ) DUMMY( : used_length )

!  Extend the length of ARRAY

         DEALLOCATE( DUMMY )
         GO TO 110
       END IF
       GO TO 10
     END IF

!  Copy the contents of ARRAY back from the buffer i/o area

       ARRAY( : used_length ) = DUMMY( : used_length )
       DEALLOCATE( DUMMY )
       new_length = length
       GO TO 200

!  Use an external unit for writing

 100   CONTINUE

!  Rewind the buffer i/o unit

     INQUIRE( UNIT = buffer, OPENED = file_open )
     IF ( file_open ) THEN
       REWIND( UNIT = buffer )
     ELSE
       OPEN( UNIT = buffer )
     END IF

!  Copy the contents of ARRAY into the buffer i/o area

     WRITE( UNIT = buffer, FMT = * ) ARRAY( : used_length )

!  Extend the length of ARRAY

     DEALLOCATE( ARRAY )

 110 CONTINUE
     ALLOCATE( ARRAY( new_length ), STAT = alloc_status )

!  If the allocation failed, reduce the new length and retry

     IF ( alloc_status /= 0 ) THEN
       new_length = min_length + ( new_length - min_length ) / 2
       IF ( new_length < min_length ) THEN
          status = 12
          RETURN
       END IF
       GO TO 110
     END IF

!  Copy the contents of ARRAY back from the buffer i/o area

     REWIND( UNIT = buffer )
     READ( UNIT = buffer, FMT = * ) ARRAY( : used_length )

!  Successful exit

 200 CONTINUE
     status = 0
     RETURN

!  end of subroutine CUTEST_extend_array_real

     END SUBROUTINE CUTEST_extend_array_real

!-  C U T E S T _ e x t e n d _ a r r a y _ i n t e g e r  S U B R O U T I N E -

     SUBROUTINE CUTEST_extend_array_integer( ARRAY, old_length, used_length,   &
                                             new_length, min_length, buffer,   &
                                             status, alloc_status )

!  -------------------------------------------------------------------------
!  extend an integer array so that its length is increaed from old_length to 
!  as close to new_length as possible while keeping existing data intact
!  -------------------------------------------------------------------------

!  History -
!   fortran 90 version released pre GALAHAD Version 1.0. February 7th 1995 as
!     EXTEND_array_integer as part of the GALAHAD module EXTEND
!   fortran 2003 version released in CUTEst, 5th November 2012

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

     INTEGER, INTENT( IN ) :: old_length, buffer
     INTEGER, INTENT( OUT ) :: status, alloc_status
     INTEGER, INTENT( INOUT ) :: used_length, min_length, new_length
     INTEGER, ALLOCATABLE, DIMENSION( : ) :: ARRAY

     INTEGER :: length
     LOGICAL :: file_open
     INTEGER, ALLOCATABLE, DIMENSION( : ) :: DUMMY

!  Make sure that the new length is larger than the old

     IF ( new_length <= old_length ) new_length = 2 * old_length

!  Ensure that the input data is consistent

     used_length = MIN( used_length, old_length )
     min_length = MAX( old_length + 1, MIN( min_length, new_length ) )

!  If possible, allocate DUMMY to hold the old values of ARRAY

     ALLOCATE( DUMMY( used_length ), STAT = alloc_status )

!  If the allocation failed, resort to using an external unit

     IF ( alloc_status /= 0 ) GO TO 100

     DUMMY( : used_length ) = ARRAY( : used_length )

!  Extend the length of ARRAY

     DEALLOCATE( ARRAY )
     length = new_length

  10 CONTINUE
     ALLOCATE( ARRAY( length ), STAT = alloc_status )

!  If the allocation failed, reduce the new length and retry

     IF ( alloc_status /= 0 ) THEN
       length = length + ( length - min_length ) / 2

!  If there is insufficient room for both ARRAY and DUMMY, use an external unit

       IF ( length < min_length ) THEN

!  Rewind the buffer i/o unit

         INQUIRE( UNIT = buffer, OPENED = file_open )
         IF ( file_open ) THEN
           REWIND( UNIT = buffer )
         ELSE
           OPEN( UNIT = buffer )
         END IF

!  Copy the contents of ARRAY into the buffer i/o area

         WRITE( UNIT = buffer, FMT = * ) DUMMY( : used_length )

!  Extend the length of ARRAY

         DEALLOCATE( DUMMY )
         GO TO 110
       END IF
       GO TO 10
     END IF

!  Copy the contents of ARRAY back from the buffer i/o area

     ARRAY( : used_length ) = DUMMY( : used_length )
     DEALLOCATE( DUMMY )
     new_length = length
     GO TO 200

!  Use an external unit for writing

 100 CONTINUE

!  Rewind the buffer i/o unit

     INQUIRE( UNIT = buffer, OPENED = file_open )
     IF ( file_open ) THEN
       REWIND( UNIT = buffer )
     ELSE
       OPEN( UNIT = buffer )
     END IF

!  Copy the contents of ARRAY into the buffer i/o area

     WRITE( UNIT = buffer, FMT = * ) ARRAY( : used_length )

!  Extend the length of ARRAY

     DEALLOCATE( ARRAY )

 110 CONTINUE
     ALLOCATE( ARRAY( new_length ), STAT = alloc_status )

!  If the allocation failed, reduce the new length and retry

     IF ( alloc_status /= 0 ) THEN
       new_length = min_length + ( new_length - min_length ) / 2
       IF ( new_length < min_length ) THEN
         status = 12
         RETURN
       END IF
       GO TO 110
     END IF

!  Copy the contents of ARRAY back from the buffer i/o area

     REWIND( UNIT = buffer )
     READ( UNIT = buffer, FMT = * ) ARRAY( : used_length )

!  Successful exit

 200 CONTINUE
     status = 0
     RETURN

!  end of subroutine CUTEST_extend_array_integer

     END SUBROUTINE CUTEST_extend_array_integer

!-*-*-*-*-*-*-*-  C U T E S T _ s y m m h  S U B R O U T I N E -*-*-*-*-*-*-*-*-

     SUBROUTINE CUTEST_symmh( maxszh, ISYMMH, ISYMMD )

!  -------------------------------------------------------------
!  Given a columnwise storage scheme of the upper triangle of a
!  symmetric matrix of order MAXSZH, compute the position of the
!  I,J-th entry of the symmetric matrix in this scheme

!  The value ISYMMH( I, J ) + 1 gives the position of the I,J-th
!  entry of the matrix in the upper triangular scheme
!  -------------------------------------------------------------

!  History -
!   fortran 77 version originally released in CUTE, September 23rd, 1991
!   fortran 90 version released pre GALAHAD Version 1.0. January 26th 1995 as
!     OTHERS_symmh as part of the GALAHAD module OTHERS
!   fortran 2003 version released in CUTEst, 5th November 2012

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

     INTEGER, INTENT( IN ) :: maxszh
     INTEGER, INTENT( OUT ), DIMENSION( maxszh, maxszh ) :: ISYMMH
     INTEGER, INTENT( OUT ), DIMENSION( maxszh ) :: ISYMMD

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

     INTEGER :: i, j, k
     
     k = 0
     DO j = 1, maxszh
       DO i = 1, j - 1
         ISYMMH( i, j ) = k ; ISYMMH( j, i ) = k ; k = k + 1
       END DO
       ISYMMD( j ) = k ; ISYMMH( j, j ) = k ; k = k + 1
     END DO
     RETURN

!  end of subroutine CUTEST_symmh

     END SUBROUTINE CUTEST_symmh

!-*-*-*-*-  C U T E S T _ g a u s s _ e l i m  S U B R O U T I N E -*-*-*-*-*-

     SUBROUTINE CUTEST_gauss_elim( m, n, IPVT, JCOL, A )

!  ------------------------------------------------------
!  Perform the first m steps of Gaussian Elimination with
!  complete pivoting on the m by n (m <= n) matrix A
!  ------------------------------------------------------

!  History -
!   fortran 77 version originally released in CUTE, September 23rd, 1991
!   fortran 90 version released pre GALAHAD Version 1.0. January 26th 1995 as
!     OTHERS_gaulss_elim as part of the GALAHAD module OTHERS
!   fortran 2003 version released in CUTEst, 5th November 2012

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

     INTEGER, INTENT( IN ) :: m, n
     INTEGER, INTENT( OUT ), DIMENSION( m ) :: IPVT
     INTEGER, INTENT( OUT ), DIMENSION( n ) :: JCOL
     REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( m * n ) :: A

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

     INTEGER :: i, j, k, l1, l2, ipivot, jpivot
     REAL ( KIND = wp ) :: apivot, atemp
     
     DO j = 1, n
       JCOL( j ) = j
     END DO

!  main loop

     DO k = 1, m

!  compute the k-th pivot

       apivot = - one
       l2 = m * ( k - 1 ) ; l1 = l2
       DO j = k, n
         DO i = k, m
           IF ( ABS( A( l1 + i ) ) > apivot ) THEN
             apivot = ABS( A( l1 + i ) )
             ipivot = i ; jpivot = j
           END IF
         END DO
         l1 = l1 + m
       END DO

!  interchange rows i and ipivot

       IPVT( k ) = ipivot
       IF ( ipivot > k ) THEN
         DO j = k, n
           atemp = A( l1 + ipivot )
           A( l2 + ipivot ) = A( l2 + k )
           A( l2 + k ) = atemp
           l2 = l2 + m
         END DO
       END IF

!  interchange columns j and jpivot

       IF ( jpivot > k ) THEN
         j = JCOL( jpivot )
         JCOL( jpivot ) = JCOL( k ) ; JCOL( k ) = j
         l1 = m * ( jpivot - 1 )
         DO i = 1, m
           atemp = A( l1 + i ) ; A( l1 + i ) = A( l2 + i )
           A( l2 + i ) = atemp
         END DO
       END IF

!  perform the elimination

       apivot = A( l2 + k )
       DO i = k + 1, m
         atemp = A( l2 + i ) / apivot
         A( l2 + i ) = atemp ; l1 = l2
         DO j = k + 1, n
           l1 = l1 + m
           A( l1 + i ) = A( l1 + i ) - atemp * A( l1 + k )
         END DO
       END DO
     END DO

     RETURN

!  end of subroutine CUTEST_gauss_elim

     END SUBROUTINE CUTEST_gauss_elim

!  end of module CUTEST

    END MODULE CUTEST
