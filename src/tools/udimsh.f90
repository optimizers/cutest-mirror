! ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE UDIMSH( data, status, nnzh )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER, INTENT( OUT ) :: status
      INTEGER :: nnzh

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

!  Local variables

      INTEGER :: nxtrw1, nxtrw2, lnxtrw, lirnh, irnh
      INTEGER :: i, ii,  ig, j,  jj, k,  l, iel, iell, inext 
      INTEGER :: newpt,  ig1,    listvs, listve, istart

!  Define the integer work space needed for ASMBI.
!  Ensure that there is sufficient space

      nnzh = 0
      newpt = data%numvar + 1
      lirnh = ( data%liwk2 - 2 * data%numvar ) / 3
      irnh = data%lsend
      lnxtrw = ( data%liwk2 - lirnh ) / 2
      nxtrw1 = irnh + lirnh
      nxtrw2 = nxtrw1 + lnxtrw
      IF ( newpt > lnxtrw .OR. lirnh <= 0 ) GO TO 900

!  NXTROW( 1, . ) gives the link list. The list for column j starts
!                 in NXTROW( 1, j ) and ends when NXTROW( 1, k ) = - 1.
!  NXTROW( 2, . ) gives the position in H of the current link.

!  Initialize the link list which points to the row numbers which
!  are used in the columns of the assembled Hessian

      DO i = 1, data%numvar
         data%IWORK( nxtrw1 + i ) = - 1
      END DO

! -------------------------------------------------------
!  Form the rank-one second order term for the IG-th group
! -------------------------------------------------------

      DO ig = 1, data%ng
        IF ( data%GXEQX( ig ) ) CYCLE
        ig1 = ig + 1
        listvs = data%IWORK( data%lstagv + ig )
        listve = data%IWORK( data%lstagv + ig1 ) - 1

!  Form the J-th column of the rank-one matrix

        DO l = listvs, listve
          j = data%IWORK( data%lsvgrp + l )
          IF ( j == 0 ) CYCLE

!  Find the entry in row i of this column

          DO k = listvs, listve
            i = data%IWORK( data%lsvgrp + k )
            IF ( i == 0 .OR. i > j ) CYCLE

!  Obtain the appropriate storage location in H for the new entry

            istart = j
  150       CONTINUE
            inext = data%IWORK( nxtrw1 + istart )
            IF ( inext == - 1 ) THEN
              IF ( newpt > lnxtrw ) GO TO 900

!  The (I,J)-th location is empty. Place the new entry in this location
!  and add another link to the list

              nnzh = nnzh + 1
              IF ( nnzh > lirnh ) GO TO 900
              data%IWORK( irnh + nnzh ) = i
              data%IWORK( nxtrw1 + istart ) = newpt
              data%IWORK( nxtrw2 + istart ) = nnzh
              data%IWORK( nxtrw1 + newpt ) = - 1
              newpt = newpt + 1

!  Continue searching the linked list for an entry in row i, column j

            ELSE
              IF ( data%IWORK( irnh + data%IWORK( nxtrw2 + istart ) ) /=i ) THEN
                istart = inext
                GO TO 150
              END IF
            END IF
          END DO
        END DO
      END DO

! --------------------------------------------------------
!  Add on the low rank first order terms for the I-th group
! --------------------------------------------------------

      DO ig = 1, data%ng
        ig1 = ig + 1

!  See if the group has any nonlinear elements

        DO iell = data%ISTADG( ig ), data%ISTADG( ig1 ) - 1
          iel = data%IELING( iell )
          listvs = data%ISTAEV( iel )
          listve = data%ISTAEV( iel + 1 ) - 1
          DO l = listvs, listve
            j = data%IELVAR( l )
            IF ( j /= 0 ) THEN

!  The IEL-th element has an internal representation.
!  Compute the J-th column of the element Hessian matrix

!  Find the entry in row i of this column

              DO k = listvs, l
                i = data%IELVAR( k )
                IF ( i /= 0 ) THEN

!  Only the upper triangle of the matrix is stored

                  IF ( i <= j ) THEN
                    ii = i
                    jj = j
                  ELSE
                    ii = j
                    jj = i
                  END IF

!  Obtain the appropriate storage location in H for the new entry

                  istart = jj
  230             CONTINUE
                  inext = data%IWORK( nxtrw1 + istart )
                  IF ( inext == - 1 ) THEN
                    IF ( newpt > lnxtrw ) GO TO 900

!  The (I,J)-th location is empty. Place the new entry in this location
!  and add another link to the list

                    nnzh = nnzh + 1
                    IF ( nnzh > lirnh ) GO TO 900
                    data%IWORK( irnh + nnzh ) = ii
                    data%IWORK( nxtrw1 + istart ) = newpt
                    data%IWORK( nxtrw2 + istart ) = nnzh
                    data%IWORK( nxtrw1 + newpt ) = - 1
                    newpt = newpt + 1
                  ELSE

!  Continue searching the linked list for an entry in row i, column j

                   IF ( data%IWORK( irnh + data%IWORK( nxtrw2 + istart ) )     &
                         == ii ) THEN
                    ELSE
                      istart = inext
                      GO TO 230
                    END IF
                  END IF
                END IF
              END DO
            END IF
          END DO
        END DO
      END DO

      status = 0
      RETURN

!  Unsuccessful returns.

  900 CONTINUE
      IF ( data%out > 0 ) WRITE( data%out,                                      &
        "( ' ** SUBROUTINE UDIMSH: Increase the size of IWK' )" )
      status = 2
      RETURN

!  end of UDIMSH.

      END
