! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE UDIMSH( data, NNZH )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: NNZH

!  Compute the space required to store the Hessian matrix of the 
!  objective function of a problem initially written in 
!  Standard Input Format (SIF).

!  NB. USETUP must have been called first

!  The upper triangle of the Hessian is stored in coordinate form,
!  i.e., the entry H(i) has row index IRNH(i)
!  for i = 1, ...., NNZH.

!  Based on the minimization subroutine data%laNCELOT/SBMIN
!  by Conn, Gould and Toint.

!  Nick Gould, for CGT productions,
!  August 1999.


! ---------------------------------------------------------------------




! ---------------------------------------------------------------------



! ---------------------------------------------------------------------


! ---------------------------------------------------------------------

!  integer variables from the GLOBAL common block.


!  integer variables from the LOCAL common block.


!  integer variables from the PRFCTS common block.


!  Local variables

      INTEGER :: NXTRW1, NXTRW2, LNXTRW, LIRNH, IRNH
      INTEGER :: I, II,  IG, J,  JJ, K,  L, IEL, IELL, INEXT 
      INTEGER :: NEWPT,  IG1,    LISTVS, LISTVE, ISTART

!  Define the integer work space needed for ASMBI.
!  Ensure that there is sufficient space

      NNZH = 0
      NEWPT = data%numvar + 1
      LIRNH = ( data%liwk2 - 2 * data%numvar ) / 3
      IRNH = data%lsend
      LNXTRW = ( data%liwk2 - LIRNH ) / 2
      NXTRW1 = IRNH + LIRNH
      NXTRW2 = NXTRW1 + LNXTRW
      IF ( NEWPT > LNXTRW .OR. LIRNH <= 0 ) GO TO 900

!  NXTROW( 1, . ) gives the link list. The list for column J starts
!                 in NXTROW( 1, J ) and ends when NXTROW( 1, K ) = - 1.
!  NXTROW( 2, . ) gives the position in H of the current link.

!  Initialize the link list which points to the row numbers which
!  are used in the columns of the assembled Hessian

      DO 20 I = 1, data%numvar
         data%IWORK( NXTRW1 + I ) = - 1
   20 CONTINUE

! -------------------------------------------------------
!  Form the rank-one second order term for the IG-th group
! -------------------------------------------------------

      DO 200 IG = 1, data%ng
         IF ( data%GXEQX( IG ) ) GO TO 200
         IG1 = IG + 1
         LISTVS = data%IWORK( data%lstagv + IG )
         LISTVE = data%IWORK( data%lstagv + IG1 ) - 1

!  Form the J-th column of the rank-one matrix

         DO 190 L = LISTVS, LISTVE
            J = data%IWORK( data%lsvgrp + L )
            IF ( J == 0 ) GO TO 190

!  Find the entry in row I of this column

            DO 180 K = LISTVS, LISTVE
               I = data%IWORK( data%lsvgrp + K )
               IF ( I == 0 .OR. I > J ) GO TO 180

!  Obtain the appropriate storage location in H for the new entry

               ISTART = J
  150          CONTINUE
               INEXT = data%IWORK( NXTRW1 + ISTART )
               IF ( INEXT == - 1 ) THEN
                  IF ( NEWPT > LNXTRW ) GO TO 900

!  The (I,J)-th location is empty. Place the new entry in this location
!  and add another link to the list

                  NNZH = NNZH + 1
                  IF ( NNZH > LIRNH ) GO TO 900
                  data%IWORK( IRNH + NNZH ) = I
                  data%IWORK( NXTRW1 + ISTART ) = NEWPT
                  data%IWORK( NXTRW2 + ISTART ) = NNZH
                  data%IWORK( NXTRW1 + NEWPT ) = - 1
                  NEWPT = NEWPT + 1
               ELSE

!  Continue searching the linked list for an entry in row I, column J

                  IF ( data%IWORK( IRNH + data%IWORK( NXTRW2 + ISTART ) )/=I ) THEN
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

      DO 300 IG = 1, data%ng
         IG1 = IG + 1

!  See if the group has any nonlinear elements

         DO 290 IELL = data%ISTADG( IG ), data%ISTADG( IG1 ) - 1
            IEL = data%IELING( IELL )
            LISTVS = data%ISTAEV( IEL )
            LISTVE = data%ISTAEV( IEL + 1 ) - 1
            DO 250 L = LISTVS, LISTVE
               J = data%IELVAR( L )
               IF ( J /= 0 ) THEN

!  The IEL-th element has an internal representation.
!  Compute the J-th column of the element Hessian matrix

!  Find the entry in row I of this column

                  DO 240 K = LISTVS, L
                     I = data%IELVAR( K )
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
                        INEXT = data%IWORK( NXTRW1 + ISTART )
                        IF ( INEXT == - 1 ) THEN
                           IF ( NEWPT > LNXTRW ) GO TO 900

!  The (I,J)-th location is empty. Place the new entry in this location
!  and add another link to the list

                           NNZH = NNZH + 1
                           IF ( NNZH > LIRNH ) GO TO 900
                           data%IWORK( IRNH + NNZH ) = II
                           data%IWORK( NXTRW1 + ISTART ) = NEWPT
                           data%IWORK( NXTRW2 + ISTART ) = NNZH
                           data%IWORK( NXTRW1 + NEWPT ) = - 1
                           NEWPT = NEWPT + 1
                        ELSE

!  Continue searching the linked list for an entry in row I, column J

                           IF ( data%IWORK( IRNH + data%IWORK( NXTRW2 + ISTART ) ) &
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

 2000 FORMAT( ' ** SUBROUTINE UDIMSH: Increase the size of IWK ' )

!  end of UDIMSH.

      END
