! THIS VERSION: CUTEST 1.0 - 20/11/2012 AT 13:30 GMT.

!-*-*-*-*-*-*-*-  C U T E S T    C F N    S U B R O U T I N E  -*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 77 version originally released in CUTE, October 1991
!   fortran 2003 version released in CUTEst, 20th November 2012

      SUBROUTINE CFN ( data, status, n, m, X, f, C )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      TYPE ( CUTEST_data_type ), INTENT( INOUT ) :: data
      INTEGER, INTENT( IN ) :: n, m
      INTEGER, INTENT( OUT ) :: status
      REAL ( KIND = wp ), INTENT( OUT ) :: f
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: X
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( m ) :: C

!  ---------------------------------------------------------------------
!  Compute the values of the objective function and general constraints
!  of a function initially written in Standard Input Format (SIF)
!  ---------------------------------------------------------------------

!  local variables

      INTEGER :: i, j, ig, ifstat, igstat
      REAL ( KIND = wp ) :: ftt

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
      END DO

!  compute the group function values

!  all group functions are trivial

      IF ( data%altriv ) THEN
        data%GVALS( : data%ng, 1 ) = data%FT( : data%ng )
        data%GVALS( : data%ng, 2 ) = 1.0_wp

!  evaluate the group function values

      ELSE
        CALL GROUP( data%GVALS, data%ng, data%FT, data%GPVALU, data%ng,        &
                    data%ITYPEG, data%ISTGP, data%ICALCF, data%ltypeg,         &
                    data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,         &
                    .FALSE., igstat )
      END IF

!  compute the objective and constraint function values.

      f = 0.0_wp
      IF ( data%numcon > 0 ) THEN
        DO ig = 1, data%ng
          IF ( data%KNDOFC( ig ) == 0 ) THEN
            IF ( data%GXEQX( ig ) ) THEN
              f = f + data%GSCALE( ig ) * data%FT( ig )
            ELSE
              f = f + data%GSCALE( ig ) * data%GVALS( ig, 1 )
            END IF
          ELSE
            IF ( data%GXEQX( ig ) ) THEN
              C( data%KNDOFC( ig ) ) = data%GSCALE( ig ) * data%FT( ig )
            ELSE
               C( data%KNDOFC( ig ) ) = data%GSCALE( ig ) * data%GVALS( ig, 1 )
            END IF
          END IF
        END DO
      ELSE

!  there are no constraints, so we need not check KNDOFC

        DO ig = 1, data%ng
          IF ( data%GXEQX( ig ) ) THEN
            f = f + data%GSCALE( ig ) * data%FT( ig )
          ELSE
            f = f + data%GSCALE( ig ) * data%GVALS( ig, 1 )
          END IF
        END DO
      END IF

!  Update the counters for the report tool.

      data%nc2of = data%nc2of + 1
      data%nc2cf = data%nc2cf + data%pnc
      status = 0
      RETURN

!  end of subroutine CFN

      END SUBROUTINE CFN
