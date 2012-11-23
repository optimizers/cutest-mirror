! THIS VERSION: CUTEST 1.0 - 23/11/2012 AT 15:00 GMT.

!-*-*-*-*-*-*-*-  C U T E S T    U S H    S U B R O U T I N E  -*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 77 version originally released in CUTE, August 1993
!   fortran 2003 version released in CUTEst, 23rd November 2012

      SUBROUTINE UBANDH( data, status, n, X, semibandwidth, H_band,            &
                         lbandh, max_semibandwidth )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      TYPE ( CUTEST_data_type ), INTENT( INOUT ) :: data
      INTEGER, INTENT( IN ) :: n, semibandwidth, lbandh
      INTEGER, INTENT( OUT ) :: status, max_semibandwidth
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: X
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( 0 : lbandh, n ) :: H_band

!  ----------------------------------------------------------------------
!  Compute the portion of the Hessian matrix of a group partially
!  separable function which lies within a band of given semi-bandwidth.
!  The diagonal and subdiagonal entries in column i are stored in
!  locations BAND( j, i ), where j = 0 for the diagonal and j > 0 for
!  the j-th subdiagonal entry, j = 1, ..., MIN( semibandwidth, n - i ) 
!  and semibandwidth is the requested semi-bandwidth. max_semibandwidth 
!  gives the full extent of the true band within semibandwidth, i.e, 
!  all entries between the bands of semi-bandwidths max_semibandwidth + 1 
!  and semibandwidth are zero
!  -----------------------------------------------------------------------

!  Based on the minimization subroutine data%laNCELOT/SBMIN
!  by Conn, Gould and Toint.

!  Nick Gould, for CGT productions.
!  August 1993.

!  Local variables

      INTEGER :: i, ig, j, nsemiw, inform, nnzh, ifstat, igstat, alloc_status
!     INTEGER :: IRNH(   1 ), ICNH(   1 )
      REAL ( KIND = wp ) :: ftt
      CHARACTER ( LEN = 80 ) :: bad_alloc = REPEAT( ' ', 80 )
      EXTERNAL :: RANGE 

!  check that there is room for the band matrix

      nsemiw = MAX( 0, MIN( semibandwidth, n - 1 ) )
      IF ( lbandh < nsemiw ) THEN
         IF ( data%out > 0 ) WRITE( data%out, 2020 )
         status = 2 ; RETURN
      END IF

!  there are non-trivial group functions

      DO i = 1, MAX( data%nel, data%ng )
        data%ICALCF( i ) = i
      END DO

!  evaluate the element function values

      CALL ELFUN( data%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,          &
                  data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,          &
                  data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,           &
                  data%lelvar, data%lntvar, data%lstadh, data%lstep,           &
                  data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,          &
                  1, ifstat )
      IF ( ifstat /= 0 ) GO TO 930

!  evaluate the element function values

      CALL ELFUN( data%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,          &
                  data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,          &
                  data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,           &
                  data%lelvar, data%lntvar, data%lstadh, data%lstep,           &
                  data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,          &
                  3, ifstat )
      IF ( ifstat /= 0 ) GO TO 930

!  compute the group argument values ft

      DO ig = 1, data%ng
        ftt = - data%B( ig )

!  include the contribution from the linear element

        DO j = data%ISTADA( ig ), data%ISTADA( ig + 1 ) - 1
          ftt = ftt + data%A( j ) * X( data%ICNA( j ) )
        END DO

!  include the contributions from the nonlinear elements

        DO j = data%ISTADG( ig ), data%ISTADG( ig + 1 ) - 1
          ftt = ftt + data%ESCALE( j ) * data%FUVALS( data%IELING( j ) )
        END DO
        data%FT( ig ) = ftt

!  Record the derivatives of trivial groups

        IF ( data%GXEQX( ig ) ) THEN
          data%GVALS( ig, 2 ) = 1.0_wp
          data%GVALS( ig, 3 ) = 0.0_wp
        END IF
      END DO

!  evaluate the group derivative values

        IF ( .NOT. data%altriv ) THEN
          CALL GROUP( data%GVALS, data%ng, data%FT, data%GPVALU, data%ng,      &
                      data%ITYPEG, data%ISTGP, data%ICALCF, data%ltypeg,       &
                      data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,       &
                      .TRUE., igstat )
          IF ( igstat /= 0 ) GO TO 930
        END IF

!  Compute the gradient value

!!$      CALL ELGRD( n, data%ng, data%firstg, data%ICNA( 1 ), data%licna, &
!!$                   data%ISTADA( 1 ), data%lstada, data%IELING( 1 ), &
!!$                   data%leling, data%ISTADG( 1 ), data%lstadg, &
!!$                   data%ITYPEE( 1 ), data%lintre, &
!!$                   data%ISTAEV( 1 ), data%lstaev, &
!!$                   data%IELVAR( 1 ), data%lelvar, &
!!$                   data%INTVAR( 1 ), data%lntvar, &
!!$                   data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, &
!!$                   data%IWORK( data%lstajc + 1 ), data%lnstjc, &
!!$                   data%IWORK( data%lstagv + 1 ), data%lnstgv, &
!!$                   data%A( 1 ), data%la, &
!!$                   data%GVALS( : , 2 ), data%lgvals, &
!!$                   data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ), &
!!$                   data%GSCALE( 1 ), data%lgscal, &
!!$                   data%ESCALE( 1 ), data%lescal, &
!!$                   data%FUVALS( data%lgrjac + 1 ), &
!!$                   data%lngrjc, data%WRK( 1 ), data%WRK( n + 1 ), data%maxsel, &
!!$                   data%GXEQX( 1 ), data%lgxeqx, &
!!$                   data%INTREP( 1 ), data%lintre, RANGE )

      CALL CUTEST_form_gradients( n, data%ng, data%nel, data%ntotel,           &
             data%nvrels, data%nnza, data%nvargp, data%firstg, data%ICNA,      &
             data%ISTADA, data%IELING, data%ISTADG, data%ISTAEV,               &
             data%IELVAR, data%INTVAR, data%A, data%GVALS( : , 2 ),            &
             data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ),          &
             data%GSCALE, data%ESCALE, data%FUVALS( data%lgrjac + 1 ),         &
             data%GXEQX, data%INTREP, data%ISVGRP, data%ISTAGV, data%ITYPEE,   &
             data%ISTAJC, data%W_ws, data%W_el, RANGE )
      data%firstg = .FALSE.

!  Define the real work space needed for ASMBL.
!  Ensure that there is sufficient space.

!     lh2 = n * ( nsemiw + 1 )
!     IF ( data%lwk2 < lh2 + n + 3 * data%maxsel ) THEN
!        IF ( data%out > 0 ) WRITE( data%out, 2000 )
!        status = 2 ; RETURN
!     END IF

!  Allocate space to hold the band matrix

!        reallocate = .TRUE.
!        IF ( ALLOCATED( DIAG ) ) THEN
!          IF ( SIZE( DIAG ) < nvar ) THEN ; DEALLOCATE( DIAG )
!           ELSE ; reallocate = .FALSE.
!           END IF
!         END IF
!         IF ( reallocate ) THEN 
!           ALLOCATE( DIAG( nvar ), STAT = alloc_status )
!           IF ( alloc_status /= 0 ) THEN ; bad_alloc = 'DIAG' ; GO TO 980
!           END IF
!         END IF
!       
!         reallocate = .TRUE.
!         IF ( ALLOCATED( OFFDIA ) ) THEN
!           IF ( SIZE( OFFDIA, 1 ) /= S%nsemiw .OR.                            &
!                SIZE( OFFDIA, 2 ) < nvar ) THEN ; DEALLOCATE( OFFDIA )
!           ELSE ; reallocate = .FALSE. ; END IF
!         END IF
!         IF ( reallocate ) THEN 
!           ALLOCATE( OFFDIA( S%nsemiw, nvar ), STAT = alloc_status )
!           IF ( alloc_status /= 0 ) THEN ; bad_alloc = 'OFFDIA'; GO TO 980
!           END IF
!         END IF

!  assemble the Hessian; use every variable

      DO i = 1, n
        data%IVAR( i ) = i
      END DO

      CALL CUTEST_assemble_hessian(                                            &
             n, data%ng, data%nel, data%ntotel, data%nvrels, data%nnza,        &
             data%maxsel, data%nvargp, n, data%IVAR, data%ISTADH,              &
             data%ICNA, data%ISTADA, data%INTVAR, data%IELVAR, data%IELING,    &
             data%ISTADG, data%ISTAEV, data%ISTAGV, data%ISVGRP, data%A,       &
             data%FUVALS, data%lnguvl, data%FUVALS, data%lnhuvl,               &
             data%GVALS( : , 2 ), data%GVALS( :  , 3 ), data%GSCALE,           &
             data%ESCALE, data%GXEQX, data%ITYPEE, data%INTREP, RANGE,         &
             0, data%out, data%out, data%io_buffer, .TRUE.,                    &
             .FALSE., .FALSE., nsemiw, status, alloc_status, bad_alloc,        &
             data%assemble_data, data%lirnh, data%ljcnh, data%lh,              &
             data%H_row, data%H_col, data%H_val,                               &
             data%LINK_col, data%POS_in_H, data%llink, data%lpos,              &
             data%IW_asmbl, data%W_ws, data%W_el, data%W_in, data%H_el,        &
             data%H_in, data%skipg, maxsbw = max_semibandwidth,                &
             DIAG = H_band( 0, : n ), OFFDIA = H_band( 1 : nsemiw, n  ) )
!            DIAG = data%DIAG, OFFDIA = data%OFFDIA )

!  check for errors in the assembly

      IF ( status > 0 ) RETURN

!  place the diagonal entries in their correct positions

!      H_band( 0, : n ) = data%DIAG( : n )
!      IF ( max_semibandwidth > 0 ) H_band( 1 : max_semibandwidth, n )         &
!        = data%OFFDIA( 1 : max_semibandwidth, n )
!      H_band( max_semibandwidth + 1, lbandh, n ) = 0.0_wp

!      CALL ASMBL( n, data%ng, data%maxsel, nsemiw, lh2, 1, nnzh, &
!                   n, data%IVAR( 1), data%ISTADH( 1 ), data%lstadh, &
!                   data%ICNA( 1 ), data%licna, &
!                   data%ISTADA( 1 ), data%lstada, &
!                   data%INTVAR( 1 ), data%lntvar, &
!                   data%IELVAR( 1 ), data%lelvar, &
!                   data%IELING( 1 ), data%leling, &
!                   data%ISTADG( 1 ), data%lstadg, &
!                   data%ISTAEV( 1 ), data%lstaev, &
!                   data%IWORK( data%lstagv + 1 ), data%lnstgv, &
!                   data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, &
!                   irnh, ICNH, data%IWORK( data%lsend + 1 ), 1, &
!                   data%IWORK( data%liwk2 - n + 1 ), n, &
!                   data%A( 1 ), data%la, &
!                   data%FUVALS, data%lnguvl, data%FUVALS, data%lnhuvl, &
!                   data%GVALS( : , 2 ), data%GVALS( : , 3 ), &
!                   data%GSCALE( 1 ), data%ESCALE( 1 ), data%lescal, &
!                   data%WRK( 1 ), data%WRK( lh2 + 1 ), &
!                   data%lwk2 - lh2, data%GXEQX( 1 ), &
!                   data%lgxeqx, data%INTREP( 1 ), data%lintre, &
!                   data%ITYPEE( 1 ), data%lintre, &
!                   RANGE, 1, data%out, .TRUE., i, inform, .FALSE., &
!                   .FALSE. )

!  Now, place the subdiagonal entries in their correct positions.

!     lh2 = n
!     DO j = 1, n
!       DO i = 1, MIN( nsemiw, n - j )
!         BANDH( i, j ) = data%WRK( lh2 + i )
!       END DO
!       lh2 = lh2 + nsemiw
!     END DO
      status = 0
      RETURN

!  unsuccessful returns

  930 CONTINUE
      IF ( data%out > 0 ) WRITE( data%out,                                     &
        "( ' ** SUBROUTINE UBANDH: error flag raised during SIF evaluation' )" )
      status = 3
      RETURN

!  non-executable statements

 2000 FORMAT( ' ** SUBROUTINE UBANDH: Increase the size of WK' )
 2020 FORMAT( ' ** SUBROUTINE UBANDH: array dimension lbandh should be',      &
              ' larger than semibandwidth' )

!  end of subroutine UBANDH

      END SUBROUTINE UBANDH

