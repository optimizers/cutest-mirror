! THIS VERSION: CUTEST 1.0 - 20/11/2012 AT 17:00 GMT.

!-*-*-*-*-*-*-*-  C U T E S T    U G R    S U B R O U T I N E  -*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 77 version originally released in CUTE, December 1990
!   fortran 2003 version released in CUTEst, 20th November 2012

      SUBROUTINE UGR( data, status, n, X, G )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      TYPE ( CUTEST_data_type ), INTENT( INOUT ) :: data
      INTEGER, INTENT( IN ) :: n
      INTEGER, INTENT( OUT ) :: status
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: X
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: G

!  ------------------------------------------------------------
!  compute the gradient of a group partially separable function
!  ------------------------------------------------------------

!  local variables

      INTEGER :: i, j, ig, ifstat, igstat
      REAL ( KIND = wp ) :: ftt
      EXTERNAL :: RANGE 

!  there are non-trivial group functions.

      DO i = 1, MAX( data%nel, data%ng )
        data%ICALCF( i ) = i
      END DO

!  evaluate the element function values.

      CALL ELFUN( data%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,          &
                  data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,          &
                  data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,           &
                  data%lelvar, data%lntvar, data%lstadh, data%lstep,           &
                  data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,          &
                  1, ifstat )
      IF ( ifstat /= 0 ) GO TO 930

! evaluate the element function derivatives

      CALL ELFUN( data%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,          &
                  data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,          &
                  data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,           &
                  data%lelvar, data%lntvar, data%lstadh, data%lstep,           &
                  data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,          &
                  2, ifstat )
      IF ( ifstat /= 0 ) GO TO 930

!  compute the group argument values ft

      DO ig = 1, data%ng
        ftt = - data%B( ig )

!  include the contribution from the linear element

        DO j = data%ISTADA( ig ), data%ISTADA( ig + 1 ) - 1
          ftt = ftt + data%A( j ) * X( data%ICNA( j ) )
        END DO

!  include the contributions from the nonlinear elements.

        DO j = data%ISTADG( ig ), data%ISTADG( ig + 1 ) - 1
          ftt = ftt + data%ESCALE( j ) * data%FUVALS( data%IELING( j ) )
        END DO
        data%FT( ig ) = ftt

!  record the derivatives of trivial groups

        IF ( data%GXEQX( ig ) ) data%GVALS( ig, 2 ) = 1.0_wp
      END DO

!  evaluate the group derivative values.

      IF ( .NOT. data%altriv ) THEN
        CALL GROUP( data%GVALS, data%ng, data%FT, data%GPVALU, data%ng,        &
                    data%ITYPEG, data%ISTGP, data%ICALCF, data%ltypeg,         &
                    data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,         &
                     .TRUE., igstat )
        IF ( igstat /= 0 ) GO TO 930
      END IF

!  Compute the gradient value.

!     CALL ELGRD( n, data%ng, data%firstg, data%ICNA, data%licna, &
!                  data%ISTADA, data%lstada, data%IELING, &
!                  data%leling, data%ISTADG, data%lstadg, &
!                  data%ITYPEE, data%lintre, &
!                  data%ISTAEV, data%lstaev, &
!                  data%IELVAR, data%lelvar, &
!                  data%INTVAR, data%lntvar, &
!                  data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, &
!                  data%IWORK( data%lstajc + 1 ), data%lnstjc, &
!                  data%IWORK( data%lstagv + 1 ), data%lnstgv, &
!                  data%A, data%la, &
!                  data%GVALS( : , 2 ), data%lgvals, &
!                  data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ), &
!                  data%GSCALE, data%lgscal, &
!                  data%ESCALE, data%lescal, &
!                  data%FUVALS( data%lgrjac + 1 ), &
!                  data%lngrjc, data%W1, data%W2, data%maxsel, &
!                  data%GXEQX, data%lgxeqx, &
!                  data%INTREP, data%lintre, RANGE )

      CALL CUTEST_form_gradients( n, data%ng, data%nel, data%ntotel,           &
             data%nvrels, data%nnza, data%nvargp, data%firstg, data%ICNA,      &
             data%ISTADA, data%IELING, data%ISTADG, data%ISTAEV,               &
             data%IELVAR, data%INTVAR, data%A, data%GVALS( : , 2 ),            &
             data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ),          &
             data%GSCALE, data%ESCALE, data%FUVALS( data%lgrjac + 1 ),         &
             data%GXEQX, data%INTREP, data%ISVGRP, data%ISTAGV, data%ITYPEE,   &
             data%ISTAJC, data%W_ws, data%W_el, RANGE )

      data%firstg = .FALSE.

!  Store the gradient value.

      DO i = 1, n
        G( i ) = data%FUVALS( data%lggfx + i )
      END DO

!  Update the counters for the report tool.

      data%nc2og = data%nc2og + 1
      status = 0
      RETURN

!  unsuccessful returns

  930 CONTINUE
      IF ( data%out > 0 ) WRITE( data%out,                                     &
        "( ' ** SUBROUTINE UGR: error flag raised during SIF evaluation' )" )
      status = 3
      RETURN

!  end of subroutine UGR

      END SUBROUTINE UGR
