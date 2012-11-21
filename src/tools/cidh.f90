! ( Last modified on 10 Sepc 2004 at 16:45:38 )
!  Correction: 10/Sep/2004: undeclared integer variables declared
      SUBROUTINE CIDH ( data, status, n, X, iprob, lh1, H )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n, iprob, lh1
      INTEGER, INTENT( OUT ) :: status
      REAL ( KIND = wp ) :: X ( n ), H( lh1, n )

!  Compute the Hessian matrix of a specified problem function 
! (IPROB = 0 is the objective function, while iprob > 0 is the 
!   IPROB-th constraint) of a problem initially written in 
!  Standard Input Format (SIF).

!  H is a two-dimensional array which gives the value of the
!    Hessian matrix of the specified problem function evaluated at
!    X. The i,j-th component of the array will contain
!    the derivative with respect to variables X(i) and X(j).

!  Based on the minimization subroutine data%laNCELOT/SBMIN
!  by Conn, Gould and Toint.

!  Nick Gould, for CGT productions,
!  August 1998.

!  Local variables

      INTEGER :: i, j, ncalcf, ncalcg, ifstat, igstat
      INTEGER :: ig, lih, liwkh, lnxtrw, linxtr, inform
      INTEGER :: k, lh, lirnh, ljcnh, nnzh, lwkh
      REAL ( KIND = wp ) :: ftt
      EXTERNAL :: RANGE 

!  Check input parameters.

      IF ( iprob < 0 ) THEN
        WRITE( data%out, 2020 )
        status = 2 ; RETURN
      END IF

!  Find group index ig of constraint IPROB.

      IF ( iprob > 0 ) THEN
        ig = 0
        DO i = 1, data%ng
          IF ( data%KNDOFC( i ) == iprob ) THEN
            ig = i
            EXIT
          END IF
        END DO
        IF ( ig == 0 ) THEN
          WRITE( data%out, 2020 )
          status = 2 ; RETURN
        END IF
      END IF

!  Find which elements are involved in the required problem function.
!  Initialize the list to zero.

      DO i = 1, data%nel
        data%ICALCF( i ) = 0
      END DO

!  If group ig is involved, record its elements.

      DO ig = 1, data%ng
        IF ( data%KNDOFC( ig ) == iprob ) THEN
          DO j = data%ISTADG( ig ), data%ISTADG( ig + 1 ) - 1
            data%ICALCF( data%IELING( j ) ) = 1
          END DO
        END IF
      END DO

!  Only compute the elements which are involved in the required function.

      ncalcf = 0
      DO i = 1, data%nel
        IF ( data%ICALCF( i ) == 1 ) THEN
          ncalcf = ncalcf + 1
          data%ICALCF( ncalcf ) = i
        ELSE

!  If this is the first ever evaluation, initialize data%FUVALS.

          IF ( data%firstg ) THEN
            data%FUVALS( i ) = 0.0_wp
            DO j = data%INTVAR( i ), data%INTVAR( i + 1 ) - 1
               data%FUVALS( j ) = 0.0_wp
            END DO
            DO j = data%ISTADH( i ), data%ISTADH( i + 1 ) - 1
               data%FUVALS( j ) = 0.0_wp
            END DO
          END IF
        END IF
      END DO

!  evaluate the element function values.

      CALL ELFUN( data%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,          &
                  data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,          &
                  data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,           &
                  data%lelvar, data%lntvar, data%lstadh, data%lstep,           &
                  data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,          &
                  1, ifstat )

!  evaluate the element function derivatives.

      CALL ELFUN( data%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,          &
                  data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,          &
                  data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,           &
                  data%lelvar, data%lntvar, data%lstadh, data%lstep,           &
                  data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,          &
                  3, ifstat )

!  Compute the list of groups involved in the required problem function.

      ncalcg = 0
      DO ig = 1, data%ng
        IF ( data%KNDOFC( ig ) == iprob ) THEN
          ncalcg = ncalcg + 1
          data%ICALCF( ncalcg ) = ig

!  compute the group argument values ft.

          ftt = - data%B( ig )

!  include the contribution from the linear element.

          DO j = data%ISTADA( ig ), data%ISTADA( ig + 1 ) - 1
            ftt = ftt + data%A( j ) * X( data%ICNA( j ) )
          END DO

!  include the contributions from the nonlinear elements.

          DO j = data%ISTADG( ig ), data%ISTADG( ig + 1 ) - 1
            ftt = ftt + data%ESCALE( j ) * data%FUVALS( data%IELING( j ) )
          END DO
          data%FT( ig ) = ftt

!  Record the derivatives of trivial groups.

          IF ( data%GXEQX( ig ) ) THEN
            data%GVALS( ig, 2 ) = 1.0_wp
            data%GVALS( ig, 3 ) = 0.0_wp
          END IF
        ELSE

!  If this is the first ever evaluation, initialize GVALS.

          IF ( data%firstg ) THEN
            data%GVALS( ig, 2 ) = 1.0_wp
            data%GVALS( ig, 3 ) = 0.0_wp
          END IF
        END IF
      END DO

!  Evaluate the group derivative values.

      IF ( .NOT. data%altriv )                                                 &
        CALL GROUP( data%GVALS, data%ng, data%FT, data%GPVALU, data%ng,        &
                    data%ITYPEG, data%ISTGP, data%ICALCF, data%ltypeg,         &
                    data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,         &
                    .TRUE., igstat )

!  Define the real work space needed for ELGRD.
!  Ensure that there is sufficient space.

      IF ( data%lwk2 < data%ng ) THEN
         IF ( data%out > 0 ) WRITE( data%out, 2000 )
         status = 2 ; RETURN
      END IF
      IF ( data%numcon > 0 ) THEN

!  Change the group weightings to include the contributions from
!  the Lagrange multipliers.

         DO ig = 1, data%ng
           i = data%KNDOFC( ig )
           IF ( i == iprob ) THEN
             data%WRK( ig ) = data%GSCALE( ig )
           ELSE
             data%WRK( ig ) = 0.0_wp
           END IF
         END DO

!  Compute the gradient value.

      CALL ELGRD( n, data%ng, data%firstg, data%ICNA( 1 ), data%licna, &
                      data%ISTADA( 1 ), data%lstada, data%IELING( 1 ), &
                      data%leling, data%ISTADG( 1 ), data%lstadg, &
                      data%ITYPEE( 1 ), data%lintre, &
                      data%ISTAEV( 1 ), data%lstaev, data%IELVAR( 1 ), &
                      data%lelvar, data%INTVAR( 1 ), data%lntvar, &
                      data%IWORK( data%lsvgrp + 1 ), &
                      data%lnvgrp, data%IWORK( data%lstajc + 1 ), data%lnstjc, &
                      data%IWORK( data%lstagv + 1 ), data%lnstgv, &
                      data%A( 1 ), data%la, &
                      data%GVALS(  : , 2 ), data%lgvals, &
                      data%FUVALS, data%lnguvl, &
                      data%FUVALS( data%lggfx + 1 ), &
                      data%WRK( 1 ), data%ng, &
                      data%ESCALE( 1 ), data%lescal, &
                      data%FUVALS( data%lgrjac + 1 ), &
                      data%lngrjc, data%WRK( 1 ), data%WRK( n + 1 ), &
                      data%maxsel, &
                      data%GXEQX( 1 ), data%lgxeqx, &
                      data%INTREP( 1 ), data%lintre, RANGE )
      ELSE

!  Compute the gradient value.

      CALL ELGRD( n, data%ng, data%firstg, data%ICNA( 1 ), data%licna, &
                      data%ISTADA( 1 ), data%lstada, data%IELING( 1 ), &
                      data%leling, data%ISTADG( 1 ), data%lstadg, &
                      data%ITYPEE( 1 ), data%lintre, &
                      data%ISTAEV( 1 ), data%lstaev, data%IELVAR( 1 ), &
                      data%lelvar, data%INTVAR( 1 ), data%lntvar, &
                      data%IWORK( data%lsvgrp + 1 ), &
                      data%lnvgrp, data%IWORK( data%lstajc + 1 ), &
                      data%lnstjc, &
                      data%IWORK( data%lstagv + 1 ), data%lnstgv, &
                      data%A( 1 ), data%la, &
                      data%GVALS( : , 2 ), data%lgvals, &
                      data%FUVALS, data%lnguvl, &
                      data%FUVALS( data%lggfx + 1 ), &
                      data%GSCALE( 1 ), data%lgscal, &
                      data%ESCALE( 1 ), data%lescal, &
                      data%FUVALS( data%lgrjac + 1 ), &
                      data%lngrjc, data%WRK( 1 ), data%WRK( n + 1 ), &
                      data%maxsel, &
                      data%GXEQX( 1 ), data%lgxeqx, &
                      data%INTREP( 1 ), data%lintre, RANGE )
      END IF
      data%firstg = .FALSE.

!  Define the real work space needed for ASMBL.
!  Ensure that there is sufficient space.

      IF ( data%numcon > 0 ) THEN
         lwkh = data%lwk2 - n - 3 * data%maxsel - data%ng
      ELSE
         lwkh = data%lwk2 - n - 3 * data%maxsel
      END IF
      IF ( lwkh <= 0 ) THEN
         IF ( data%out > 0 ) WRITE( data%out, 2000 )
         status = 2 ; RETURN
      END IF

!  Define the integer work space needed for ASMBL.
!  Ensure that there is sufficient space.

      liwkh = data%liwk2 - n
      lh = MIN( lwkh, ( liwkh - 3 * n ) / 4 )
      linxtr = lh + n
      IF ( lh <= 0 ) THEN
         IF ( data%out > 0 ) WRITE( data%out, 2010 )
         status = 2 ; RETURN
      END IF

!  Set starting addresses for partitions of the integer workspace.

      lih = lh
      lirnh = 0
      ljcnh = lirnh + lih
      lnxtrw = ljcnh + lih
      DO 90 i = 1, n
         data%IVAR( i ) = i
   90 CONTINUE

!  Assemble the Hessian.

      IF ( data%numcon > 0 ) THEN
      CALL ASMBL( n, data%ng, data%maxsel, n, lh, lih, nnzh, &
                   n, data%IVAR( 1), data%ISTADH( 1 ), data%lstadh, &
                   data%ICNA( 1 ), data%licna, &
                   data%ISTADA( 1 ), data%lstada, &
                   data%INTVAR( 1 ), data%lntvar, &
                   data%IELVAR( 1 ), data%lelvar, &
                   data%IELING( 1 ), data%leling, &
                   data%ISTADG( 1 ), data%lstadg, &
                   data%ISTAEV( 1 ), data%lstaev, &
                   data%IWORK( data%lstagv + 1 ), data%lnstgv, &
                   data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, &
                   data%IWORK( data%lsend + lirnh + 1 ), &
                   data%IWORK( data%lsend + ljcnh + 1 ), &
                   data%IWORK( data%lsend + lnxtrw + 1 ), linxtr, &
                   data%IWORK( data%lsend + liwkh + 1 ), n, &
                   data%A( 1 ), data%la, data%FUVALS, data%lnguvl, &
                   data%FUVALS, data%lnhuvl, &
                   data%GVALS( : , 2 ), data%GVALS( : , 3 ), &
                   data%WRK( 1 ), data%ESCALE( 1 ), data%lescal, &
                   data%WRK( data%ng + 1 ), data%WRK( lwkh + data%ng + 1 ), &
                   data%lwk2 - lwkh, data%GXEQX( 1 ), &
                   data%lgxeqx, data%INTREP( 1 ), data%lintre, &
                   data%ITYPEE( 1 ), data%lintre, &
                   RANGE, 1, data%out, .FALSE., i, inform, .FALSE., &
                   .FALSE. )
      ELSE
      CALL ASMBL( n, data%ng, data%maxsel, n, lh, lih, nnzh, &
                   n, data%IVAR( 1), data%ISTADH( 1 ), data%lstadh, &
                   data%ICNA( 1 ), data%licna, &
                   data%ISTADA( 1 ), data%lstada, &
                   data%INTVAR( 1 ), data%lntvar, &
                   data%IELVAR( 1 ), data%lelvar, &
                   data%IELING( 1 ), data%leling, &
                   data%ISTADG( 1 ), data%lstadg, &
                   data%ISTAEV( 1 ), data%lstaev, &
                   data%IWORK( data%lstagv + 1 ), &
                   data%lnstgv, data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, &
                   data%IWORK( data%lsend + lirnh + 1 ), &
                   data%IWORK( data%lsend + ljcnh + 1 ), &
                   data%IWORK( data%lsend + lnxtrw + 1 ), linxtr, &
                   data%IWORK( data%lsend + liwkh + 1 ), n, &
                   data%A( 1 ), data%la, &
                   data%FUVALS, data%lnguvl, data%FUVALS, data%lnhuvl, &
                   data%GVALS( : , 2 ), data%GVALS( : , 3 ), &
                   data%GSCALE( 1 ), data%ESCALE( 1 ), data%lescal, &
                   data%WRK( 1 ), data%WRK( lwkh + 1 ), &
                   data%lwk2 - lwkh, data%GXEQX( 1 ), &
                   data%lgxeqx, data%INTREP( 1 ), data%lintre, &
                   data%ITYPEE( 1 ), data%lintre, &
                   RANGE, 1, data%out, .FALSE., i, inform, .FALSE., &
                   .FALSE. )
      END IF

!  Check that there is sufficient integer workspace.

      IF ( inform > 0 ) THEN
         IF ( data%out > 0 ) WRITE( data%out, 2010 )
         status = 2 ; RETURN
      END IF

!  Initialize the dense matrix.

      H( : n, : n ) = 0.0_wp

!  Transfer the matrix from co-ordinate to dense storage and
!  symmetrize the martix.

      IF ( data%numcon > 0 ) THEN
        DO k = 1, nnzh
          i = data%IWORK( data%lsend + lirnh + k )
          j = data%IWORK( data%lsend + ljcnh + k )
          H( i, j ) = data%WRK( data%ng + k )
          H( j, i ) = data%WRK( data%ng + k )
        END DO
      ELSE
        DO k = 1, nnzh
          i = data%IWORK( data%lsend + lirnh + k )
          j = data%IWORK( data%lsend + ljcnh + k )
          H( i, j ) = data%WRK( k )
          H( j, i ) = data%WRK( k )
        END DO
      END IF

!  Update the counters for the report tool.

      IF ( iprob == 0 ) THEN
        data%nc2oh = data%nc2oh + 1
      ELSE 
        data%nc2ch = data%nc2ch + 1
      ENDIF 
      status = 0
      RETURN

! Non-executable statements.

 2000 FORMAT( ' ** SUBROUTINE CIDH: Increase the size of WK ' )
 2010 FORMAT( ' ** SUBROUTINE CIDH: Increase the size of IWK ' )
 2020 FORMAT( ' ** SUBROUTINE CIDH: invalid problem index iprob ' )

!  end of CIDH.

      END
