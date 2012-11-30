! THIS VERSION: CUTEST 1.0 - 28/11/2012 AT 10:15 GMT.

!-*-*-*-*-*-*-  C U T E S T    C C I F S G    S U B R O U T I N E  -*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal authors: Ingrid Bongartz and Nick Gould

!   replaced obsolete subroutine CSCIFG in CUTEr, September 1994
!   fortran 2003 version released in CUTEst, 28th November 2012

      SUBROUTINE CCIFSG( data, status, n, icon, X, ci,                         &
                         nnzgci, lgci, GCI_val, GCI_var, grad )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      TYPE ( CUTEST_data_type ), INTENT( INOUT ) :: data
      INTEGER, INTENT( IN ) :: n, icon, lgci
      INTEGER, INTENT( OUT ) :: status, nnzgci
      LOGICAL, INTENT( IN ) :: grad
      INTEGER, INTENT( OUT ), DIMENSION( lgci ) :: GCI_var
      REAL ( KIND = wp ), INTENT( OUT ) :: ci
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: X
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( lgci ) :: GCI_val

!  ---------------------------------------------------------------------
!  evaluate constraint function icon and possibly its gradient,
!  for constraints initially written in Standard Input Format (SIF).
!  The constraint gradient is stored as a sparse vector in array GCI_val.
!  The j-th entry of GCI_val gives the value of the partial derivative
!  of constraint icon with respect to variable GCI_var(j).
!  The number of nonzeros in vector GCI_val is given by nnzgci.
!  (Subroutine CCIFG performs the same calculations for a dense
!   constraint gradient vector.)
!  ---------------------------------------------------------------------

!  local variables

      INTEGER :: i, j, iel, k, ig, ii, ig1, l, ll, neling
      INTEGER :: nin, nvarel, nelow, nelup, istrgv, iendgv, ifstat, igstat
      EXTERNAL :: RANGE 
      REAL ( KIND = wp ) :: ftt, gi, scalee

!  return if there are no constraints

      IF ( data%numcon == 0 ) RETURN

!  check input parameters

      IF ( icon <= 0 ) THEN
        IF ( data%out > 0 ) WRITE( data%out, 2000 )
        status = 2 ; RETURN
      END IF

!  find group index ig of constraint icon

      ig = 0
      DO i = 1, data%ng
        IF ( data%KNDOFC( i ) == icon ) THEN
          ig = i
          EXIT
        END IF
      END DO
      IF ( ig == 0 ) THEN
        IF ( data%out > 0 ) WRITE( data%out, 2000 )
        status = 2 ; RETURN
      END IF

!  determine nonlinear elements in group ig. Record their indices in ICALCF

      nelow = data%ISTADG( ig )
      nelup = data%ISTADG( ig + 1 ) - 1
      neling = nelup - nelow + 1
      j = nelow - 1
      DO i = 1, neling
        j = j + 1
        data%ICALCF( i ) = data%IELING( j )
      END DO

!  evaluate the element function values

      CALL ELFUN( data%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,          &
                  data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,          &
                  data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,           &
                  data%lelvar, data%lntvar, data%lstadh, data%lstep,           &
                  data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,          &
                  1, ifstat )
      IF ( ifstat /= 0 ) GO TO 930

!  compute the group argument value FTT. Consider only the group associated 
!  with constraint icon

      ftt = - data%B( ig )

!  include contributions from the linear element only if the variable belongs 
!  to the first n variables

      DO i = data%ISTADA( ig ), data%ISTADA( ig + 1 ) - 1
        j = data%ICNA( i )
        IF ( j <= n ) ftt = ftt + data%A( i ) * X( j )
      END DO

!  include the contributions from the nonlinear elements

      DO i = nelow, nelup
        ftt = ftt + data%ESCALE( i ) * data%FUVALS( data%IELING( i ) )
      END DO
      data%FT( ig ) = ftt

!  if ig is a trivial group, record the function value and derivative.

      IF ( data%GXEQX( ig ) ) THEN
        data%GVALS( ig, 1 ) = data%FT( ig )
        data%GVALS( ig, 2 ) = 1.0_wp

!  Otherwise, evaluate group IG.

      ELSE
        CALL GROUP( data%GVALS, data%ng, data%FT, data%GPVALU, data%ng,        &
                    data%ITYPEG, data%ISTGP, data%ICALCF, data%ltypeg,         &
                    data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,         &
                    .FALSE., igstat )
        IF ( igstat /= 0 ) GO TO 930
      END IF

!  compute the constraint function value

      IF ( data%GXEQX( ig ) ) THEN
        ci = data%GSCALE( ig ) * data%FT( ig )
      ELSE 
        ci = data%GSCALE( ig ) * data%GVALS( ig, 1 )
      END IF

!  evaluate the element function derivatives

      IF ( grad ) THEN
        CALL ELFUN( data%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,        &
                    data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,        &
                    data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,         &
                    data%lelvar, data%lntvar, data%lstadh, data%lstep,         &
                    data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,        &
                    2, ifstat )
        IF ( ifstat /= 0 ) GO TO 930

!  evaluate the group derivative values

        IF ( .NOT. data%GXEQX( ig ) ) THEN
          CALL GROUP( data%GVALS, data%ng, data%FT, data%GPVALU, data%ng,      &
                      data%ITYPEG, data%ISTGP, data%ICALCF, data%ltypeg,       &
                      data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,       &
                      .TRUE., igstat )
          IF ( igstat /= 0 ) GO TO 930
        END IF

!  compute the gradient values. Initialize the gradient vector as zero

        nnzgci = 0 
        GCI_val( : lgci ) = 0.0_wp

!  consider only group ig

        ig1 = ig + 1
        istrgv = data%ISTAGV( ig )
        iendgv = data%ISTAGV( ig1 ) - 1

!  compute the first derivative of the group

        gi = data%GSCALE( ig )
        IF ( .NOT. data%GXEQX( ig ) ) gi = gi  * data%GVALS( ig, 2 )
        data%W_ws( data%ISVGRP( istrgv : iendgv ) ) = 0.0_wp

!  the group has nonlinear elements

        IF ( nelow <= nelup ) THEN

!  loop over the group's nonlinear elements

          DO ii = nelow, nelup
            iel = data%IELING( ii )
            k = data%INTVAR( iel )
            l = data%ISTAEV( iel )
            nvarel = data%ISTAEV( iel + 1 ) - l
            scalee = data%ESCALE( ii )
            IF ( data%INTREP( iel ) ) THEN

!  the iel-th element has an internal representation

              nin = data%INTVAR( iel + 1 ) - k
              CALL RANGE( iel, .TRUE., data%FUVALS( k ), data%W_el,            &
                          nvarel, nin, data%ITYPEE( iel ), nin, nvarel )
!DIR$ IVDEP
              DO i = 1, nvarel
                j = data%IELVAR( l )
                data%W_ws( j ) = data%W_ws( j ) + scalee * data%W_el( i )
                l = l + 1
              END DO
            ELSE

!  the iel-th element has no internal representation

!DIR$ IVDEP
              DO i = 1, nvarel
                j = data%IELVAR( l )
                data%W_ws( j ) = data%W_ws( j ) + scalee * data%FUVALS( k )
                k = k + 1 ; l = l + 1
              END DO
            END IF
          END DO

!  include the contribution from the linear element

!DIR$ IVDEP
          DO k = data%ISTADA( ig ), data%ISTADA( ig1 ) - 1
               j = data%ICNA( k )
               data%W_ws( j ) = data%W_ws( j ) + data%A( k )
          END DO

!  allocate a gradient

!DIR$ IVDEP
          DO i = istrgv, iendgv
            ll = data%ISVGRP( i )

!  include contributions from the first n variables only

            IF ( ll <= n ) THEN
              nnzgci = nnzgci + 1
              GCI_val ( nnzgci ) = gi * data%W_ws( ll )
              GCI_var( nnzgci ) = ll
            END IF
          END DO 

!  the group has only linear elements

        ELSE

!  include the contribution from the linear element

!DIR$ IVDEP
          DO k = data%ISTADA( ig ),data%ISTADA( ig1 ) - 1
            j = data%ICNA( k )
            data%W_ws( j ) = data%W_ws( j ) + data%A( k )
          END DO

!  allocate a gradient

!DIR$ IVDEP
          DO i = istrgv, iendgv
            ll = data%ISVGRP( i )

!  include contributions from the first n variables only

            IF ( ll <= n ) THEN
              nnzgci = nnzgci + 1
              GCI_val ( nnzgci ) = gi * data%W_ws( ll )
              GCI_var( nnzgci ) = ll
            END IF
          END DO
        END IF
      END IF

!  update the counters for the report tool

      data%nc2cf = data%nc2cf + 1
      IF ( grad ) data%nc2cg = data%nc2cg + 1
      status = 0
      RETURN

!  unsuccessful returns

  930 CONTINUE
      IF ( data%out > 0 ) WRITE( data%out,                                     &
        "( ' ** SUBROUTINE CCIFSG: error flag raised during SIF evaluation' )" )
      status = 3
      RETURN

!  non-executable statements

 2000 FORMAT( ' ** SUBROUTINE CCIFSG: invalid constraint index icon ' )

!  end of subroutine CCIFSG

      END SUBROUTINE CCIFSG
