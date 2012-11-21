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

! evaluate the element function derivatives

      CALL ELFUN( data%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,          &
                  data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,          &
                  data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,           &
                  data%lelvar, data%lntvar, data%lstadh, data%lstep,           &
                  data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,          &
                  2, ifstat )

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

      IF ( .NOT. data%altriv )                                                 &
        CALL GROUP( data%GVALS, data%ng, data%FT, data%GPVALU, data%ng,        &
                    data%ITYPEG, data%ISTGP, data%ICALCF, data%ltypeg,         &
                    data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,         &
                     .TRUE., igstat )

!  Compute the gradient value.

      CALL ELGRD( n, data%ng, data%firstg, data%ICNA( 1 ), data%licna, &
                   data%ISTADA( 1 ), data%lstada, data%IELING( 1 ), &
                   data%leling, data%ISTADG( 1 ), data%lstadg, &
                   data%ITYPEE( 1 ), data%lintre, &
                   data%ISTAEV( 1 ), data%lstaev, &
                   data%IELVAR( 1 ), data%lelvar, &
                   data%INTVAR( 1 ), data%lntvar, &
                   data%IWORK( data%lsvgrp + 1 ), data%lnvgrp, &
                   data%IWORK( data%lstajc + 1 ), data%lnstjc, &
                   data%IWORK( data%lstagv + 1 ), data%lnstgv, &
                   data%A( 1 ), data%la, &
                   data%GVALS( : , 2 ), data%lgvals, &
                   data%FUVALS, data%lnguvl, data%FUVALS( data%lggfx + 1 ), &
                   data%GSCALE( 1 ), data%lgscal, &
                   data%ESCALE( 1 ), data%lescal, &
                   data%FUVALS( data%lgrjac + 1 ), &
                   data%lngrjc, data%WRK( 1 ), data%WRK( n + 1 ), data%maxsel, &
                   data%GXEQX( 1 ), data%lgxeqx, &
                   data%INTREP( 1 ), data%lintre, RANGE )
      data%firstg = .FALSE.

!  Store the gradient value.

      DO i = 1, n
        G( i ) = data%FUVALS( data%lggfx + i )
      END DO

!  Update the counters for the report tool.

      data%nc2og = data%nc2og + 1
      status = 0
      RETURN

!  end of subroutine UGR

      END SUBROUTINE UGR
