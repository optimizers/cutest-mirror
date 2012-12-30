! THIS VERSION: CUTEST 1.0 - 29/12/2012 AT 13:15 GMT.

!-*-*-*-*-*-*-  C U T E S T    C C I F S G    S U B R O U T I N E  -*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 2003 version released in CUTEst, 29th December 2012

      SUBROUTINE CUTEST_ccifsg( status, n, icon, X, ci,                        &
                                nnzgci, lgci, GCI_val, GCI_var, grad )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

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

      CALL CUTEST_ccifsg_threadsafe( CUTEST_data_global,                       &
                                     CUTEST_work_global( 1 ),                  &
                                     status, n, icon, X, ci, nnzgci, lgci,     &
                                     GCI_val, GCI_var, grad )
      RETURN

!  end of subroutine CUTEST_ccifsg

      END SUBROUTINE CUTEST_ccifsg

!-*-  C U T E S T    C C I F S G _ t h r e a d e d   S U B R O U T I N E  -*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 2003 version released in CUTEst, 29th December 2012

      SUBROUTINE CUTEST_ccifsg_threaded( status, n, icon, X, ci, nnzgci, lgci, &
                                         GCI_val, GCI_var, grad, thread )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      INTEGER, INTENT( IN ) :: n, icon, lgci, thread
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

      CALL CUTEST_ccifsg_threadsafe( CUTEST_data_global,                       &
                                     CUTEST_work_global( thread ),             &
                                     status, n, icon, X, ci, nnzgci, lgci,     &
                                     GCI_val, GCI_var, grad )
      RETURN

!  end of subroutine CUTEST_ccifsg_threaded

      END SUBROUTINE CUTEST_ccifsg_threaded

!-*-  C U T E S T   C C I F S G _ t h r e a d s a f e   S U B R O U T I N E  -*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal authors: Ingrid Bongartz and Nick Gould

!  History -
!   replaced obsolete subroutine CSCIFG in CUTEr, September 1994
!   fortran 2003 version released in CUTEst, 28th November 2012

      SUBROUTINE CUTEST_ccifsg_threadsafe( data, work, status, n, icon, X,     &
                                           ci, nnzgci, lgci, GCI_val,          &
                                           GCI_var, grad )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      TYPE ( CUTEST_data_type ), INTENT( IN ) :: data
      TYPE ( CUTEST_work_type ), INTENT( INOUT ) :: work
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
        work%ICALCF( i ) = data%IELING( j )
      END DO

!  evaluate the element function values

      CALL ELFUN( work%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,          &
                  data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,          &
                  data%ISTEP, work%ICALCF, data%ltypee, data%lstaev,           &
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
        ftt = ftt + data%ESCALE( i ) * work%FUVALS( data%IELING( i ) )
      END DO
      work%FT( ig ) = ftt

!  if ig is a trivial group, record the function value and derivative.

      IF ( data%GXEQX( ig ) ) THEN
        work%GVALS( ig, 1 ) = work%FT( ig )
        work%GVALS( ig, 2 ) = 1.0_wp

!  Otherwise, evaluate group IG.

      ELSE
        CALL GROUP( work%GVALS, data%ng, work%FT, data%GPVALU, data%ng,        &
                    data%ITYPEG, data%ISTGP, work%ICALCF, data%ltypeg,         &
                    data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,         &
                    .FALSE., igstat )
        IF ( igstat /= 0 ) GO TO 930
      END IF

!  compute the constraint function value

      IF ( data%GXEQX( ig ) ) THEN
        ci = data%GSCALE( ig ) * work%FT( ig )
      ELSE 
        ci = data%GSCALE( ig ) * work%GVALS( ig, 1 )
      END IF

!  evaluate the element function derivatives

      IF ( grad ) THEN
        CALL ELFUN( work%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,        &
                    data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,        &
                    data%ISTEP, work%ICALCF, data%ltypee, data%lstaev,         &
                    data%lelvar, data%lntvar, data%lstadh, data%lstep,         &
                    data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,        &
                    2, ifstat )
        IF ( ifstat /= 0 ) GO TO 930

!  evaluate the group derivative values

        IF ( .NOT. data%GXEQX( ig ) ) THEN
          CALL GROUP( work%GVALS, data%ng, work%FT, data%GPVALU, data%ng,      &
                      data%ITYPEG, data%ISTGP, work%ICALCF, data%ltypeg,       &
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
        IF ( .NOT. data%GXEQX( ig ) ) gi = gi  * work%GVALS( ig, 2 )
        work%W_ws( data%ISVGRP( istrgv : iendgv ) ) = 0.0_wp

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
              CALL RANGE( iel, .TRUE., work%FUVALS( k ), work%W_el,            &
                          nvarel, nin, data%ITYPEE( iel ), nin, nvarel )
!DIR$ IVDEP
              DO i = 1, nvarel
                j = data%IELVAR( l )
                work%W_ws( j ) = work%W_ws( j ) + scalee * work%W_el( i )
                l = l + 1
              END DO
            ELSE

!  the iel-th element has no internal representation

!DIR$ IVDEP
              DO i = 1, nvarel
                j = data%IELVAR( l )
                work%W_ws( j ) = work%W_ws( j ) + scalee * work%FUVALS( k )
                k = k + 1 ; l = l + 1
              END DO
            END IF
          END DO

!  include the contribution from the linear element

!DIR$ IVDEP
          DO k = data%ISTADA( ig ), data%ISTADA( ig1 ) - 1
               j = data%ICNA( k )
               work%W_ws( j ) = work%W_ws( j ) + data%A( k )
          END DO

!  allocate a gradient

!DIR$ IVDEP
          DO i = istrgv, iendgv
            ll = data%ISVGRP( i )

!  include contributions from the first n variables only

            IF ( ll <= n ) THEN
              nnzgci = nnzgci + 1
              GCI_val ( nnzgci ) = gi * work%W_ws( ll )
              GCI_var( nnzgci ) = ll
            END IF
          END DO 

!  the group has only linear elements

        ELSE

!  include the contribution from the linear element

!DIR$ IVDEP
          DO k = data%ISTADA( ig ),data%ISTADA( ig1 ) - 1
            j = data%ICNA( k )
            work%W_ws( j ) = work%W_ws( j ) + data%A( k )
          END DO

!  allocate a gradient

!DIR$ IVDEP
          DO i = istrgv, iendgv
            ll = data%ISVGRP( i )

!  include contributions from the first n variables only

            IF ( ll <= n ) THEN
              nnzgci = nnzgci + 1
              GCI_val ( nnzgci ) = gi * work%W_ws( ll )
              GCI_var( nnzgci ) = ll
            END IF
          END DO
        END IF
      END IF

!  update the counters for the report tool

      work%nc2cf = work%nc2cf + 1
      IF ( grad ) work%nc2cg = work%nc2cg + 1
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

!  end of subroutine CUTEST_ccifsg_threadsafe

      END SUBROUTINE CUTEST_ccifsg_threadsafe
