! THIS VERSION: CUTEST 1.0 - 21/11/2012 AT 08:30 GMT.

!-*-*-*-*-*-*-*-  C U T E S T    C C F G    S U B R O U T I N E  -*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal authors: Ingrid Bongartz and Nick Gould

!  History -
!   fortran 77 version originally released in CUTE, April 1992
!   fortran 2003 version released in CUTEst, 21st November 2012

      SUBROUTINE CCFG( data, status, n, m, X, C, jtrans,                       &
                       lcjac1, lcjac2, CJAC, grad )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      TYPE ( CUTEST_data_type ), INTENT( INOUT ) :: data
      INTEGER, INTENT( IN ) :: n, m, lcjac1, lcjac2
      INTEGER, INTENT( OUT ) :: status
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: X
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( m ) :: C
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( lcjac1, lcjac2 ) :: CJAC
      LOGICAL, INTENT( IN ) :: jtrans, grad

! ------------------------------------------------------------------------
!  Compute the values of the constraint functions and their gradients
!  for constraints initially written in Standard Input Format (SIF).  
!  The Jacobian must be stored in a dense format.
! (Subroutine CSCFG performs the same calculations for a sparse Jacobian.)

!  CJAC  is a two-dimensional array of dimension ( lcjac1, lcjac2 )
!        which gives the value of the Jacobian matrix of the
!        constraint functions, or its transpose, evaluated at X.
!        If JTRANS is .TRUE., the i,j-th component of the array
!        will contain the i-th derivative of the j-th constraint
!        function. Otherwise, if JTRANS is .FALSE., the i,j-th
!        component of the array will contain the j-th derivative
!        of the i-th constraint function.
! ------------------------------------------------------------------------

!  local variables

      INTEGER :: i, j, iel, k, ig, ii, ig1, l, ll, icon, nin, nvarel, nelow
      INTEGER :: icnt, ifstat, igstat, nelup, istrgv, iendgv
      EXTERNAL :: RANGE 
      REAL ( KIND = wp ) :: ftt, gi, scalee

      IF ( data%numcon == 0 ) RETURN

!  check input parameters

      IF ( grad ) THEN
        IF ( jtrans ) THEN
          IF ( lcjac1 < n .OR. lcjac2 < m ) THEN
            IF ( lcjac1 < n .AND. data%out > 0 ) WRITE( data%out, 2000 )
            IF ( lcjac2 < m .AND. data%out > 0 ) WRITE( data%out, 2010 )
            status = 2 ; RETURN
          END IF
        ELSE
          IF ( lcjac1 < m .OR. lcjac2 < n ) THEN
            IF ( lcjac1 < m .AND. data%out > 0 ) WRITE( data%out, 2000 )
            IF ( lcjac2 < n .AND. data%out > 0 ) WRITE( data%out, 2010 )
            status = 2 ; RETURN
          END IF
        END IF
      END IF

!  identify which elements are included in constraints. Use logical work 
!  vector to keep track of elements already included

      data%LOGIC( 1 : data%nel ) = .FALSE.

!  Now identify elements in first m constraint groups

      icnt = 0
      DO ig = 1, data%ng
        icon = data%KNDOFC( ig )
        IF ( icon > 0 .AND. icon <= m ) THEN
          nelow = data%ISTADG( ig )
          nelup = data%ISTADG( ig + 1 ) - 1
          DO ii = nelow, nelup
            iel = data%IELING( ii )
            IF ( .NOT. data%LOGIC( iel ) ) THEN
              data%LOGIC( iel ) = .TRUE.
              icnt = icnt + 1
              data%ICALCF( icnt ) = iel
            END IF
          END DO
        END IF
      END DO

!  evaluate the element function values

      CALL ELFUN( data%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,          &
                  data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,          &
                  data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,           &
                  data%lelvar, data%lntvar, data%lstadh, data%lstep,           &
                  data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,          &
                  1, ifstat )
      IF ( ifstat /= 0 ) GO TO 930

!  evaluate the element function derivatives

      IF ( GRAD )                                                              &
        CALL ELFUN( data%FUVALS, X, data%EPVALU, data%nel, data%ITYPEE,        &
                    data%ISTAEV, data%IELVAR, data%INTVAR, data%ISTADH,        &
                    data%ISTEP, data%ICALCF, data%ltypee, data%lstaev,         &
                    data%lelvar, data%lntvar, data%lstadh, data%lstep,         &
                    data%lcalcf, data%lfuval, data%lvscal, data%lepvlu,        &
                    2, ifstat )
      IF ( ifstat /= 0 ) GO TO 930

!  compute the group argument values ft

      DO ig = 1, data%ng
         ftt = 0.0_wp

!  consider only those groups in the constraints

         icon = data%KNDOFC( ig )
         IF ( icon > 0 .AND. icon <= m ) THEN
           ftt = - data%B( ig )

!  include the contribution from the linear element only if the variable 
!  belongs to the first n variables

           DO i = data%ISTADA( ig ), data%ISTADA( ig + 1 ) - 1
             j = data%ICNA( i )
             IF ( j <= n ) ftt = ftt + data%A( i ) * X( j )
           END DO

!  include the contributions from the nonlinear elements

           DO i = data%ISTADG( ig ), data%ISTADG( ig + 1 ) - 1
             ftt = ftt + data%ESCALE( i ) * data%FUVALS( data%IELING( i ) )
           END DO

!  record the derivatives of trivial groups

           IF ( data%GXEQX( ig ) ) data%GVALS( ig, 2 ) = 1.0_wp
         END IF
         data%FT( ig ) = ftt
       END DO

!  compute the group function values

!  all group functions are trivial

      IF ( data%altriv ) THEN
        data%GVALS( : data%ng, 1 ) = data%FT( : data%ng )
        data%GVALS( : data%ng, 2 ) = 1.0_wp

!  evaluate the group function values. Evaluate groups belonging to the first 
!  m constraints only

      ELSE
        icnt = 0
        DO ig = 1, data%ng
          icon = data%KNDOFC( ig )
          IF ( icon > 0 .AND. icon <= m ) THEN
            icnt = icnt + 1
            data%ICALCF( icnt ) = ig
          END IF 
        END DO
        CALL GROUP( data%GVALS, data%ng, data%FT, data%GPVALU, data%ng,        &
                    data%ITYPEG, data%ISTGP, data%ICALCF, data%ltypeg,         &
                    data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,         &
                    .FALSE., igstat )
        IF ( igstat /= 0 ) GO TO 930
      END IF

!  compute the constraint function values

      DO ig = 1, data%ng
        i = data%KNDOFC( ig )
        IF ( i > 0 .AND. i <= m ) THEN
          IF ( data%GXEQX( ig ) ) THEN
            C( i ) = data%GSCALE( ig ) * data%FT( ig )
          ELSE
            C( i ) = data%GSCALE( ig ) * data%GVALS( ig, 1 )
          END IF
        END IF
      END DO

!  increment the constraint function evaluation counter

      data%nc2cf = data%nc2cf + data%pnc

!  increment the constraint gradient evaluation counter

      IF ( GRAD ) THEN
        data%nc2cg = data%nc2cg + data%pnc

!  evaluate the group derivative values

        IF ( .NOT. data%altriv ) THEN
          CALL GROUP( data%GVALS, data%ng, data%FT, data%GPVALU, data%ng,      &
                      data%ITYPEG, data%ISTGP, data%ICALCF, data%ltypeg,       &
                      data%lstgp, data%lcalcf, data%lcalcg, data%lgpvlu,       &
                      .TRUE., igstat )
          IF ( igstat /= 0 ) GO TO 930
        END IF

!  compute the gradient values.  Initialize the Jacobian as zero

        IF ( JTRANS ) THEN
           CJAC( 1 : n, 1 : m ) = 0.0_wp
        ELSE
           CJAC( 1 : m, 1 : n ) = 0.0_wp
        END IF

!  consider the ig-th group

        DO ig = 1, data%ng
          icon = data%KNDOFC( ig )

!  consider only those groups in the first m constraints

          IF ( icon == 0 .OR. icon > m ) CYCLE
          ig1 = ig + 1
          istrgv = data%ISTAGV( ig )
          iendgv = data%ISTAGV( ig1 ) - 1
          nelow = data%ISTADG( ig ) ; nelup = data%ISTADG( ig1 ) - 1

!  compute the first derivative of the group

          gi = data%GSCALE( ig )
          IF ( .NOT. data%GXEQX( ig ) ) gi = gi  * data%GVALS( ig, 2 )

!  the group has nonlinear elements

          IF ( nelow <= nelup ) THEN
            data%W_ws( data%ISVGRP( istrgv : iendgv ) ) = 0.0_wp

!  loop over the group's nonlinear elements

            DO ii = nelow, nelup
              iel = data%IELING( ii )
              k = data%INTVAR( iel ) ; l = data%ISTAEV( iel )
              nvarel = data%ISTAEV( iel + 1 ) - l
              scalee = data%ESCALE( ii )

!  the iel-th element has an internal representation

              IF ( data%INTREP( iel ) ) THEN
                nin = data%INTVAR( iel + 1 ) - k
                CALL RANGE( iel, .TRUE., data%FUVALS( k ), data%W_el,          &
                            nvarel, nin, data%ITYPEE( iel ), nin, nvarel )
!DIR$ IVDEP
                DO i = 1, nvarel
                  j = data%IELVAR( l )
                  data%W_ws( j ) = data%W_ws( j ) + scalee * data%W_el( i )
                  l = l + 1
                END DO

!  the iel-th element has no internal representation

              ELSE
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
                IF ( jtrans ) THEN
                  CJAC( ll, icon ) = gi * data%W_ws( ll )
                ELSE
                  CJAC( icon, ll ) = gi * data%W_ws( ll )
                END IF
              END IF
            END DO

!  the group has only linear elements

          ELSE

!  allocate a gradient

!DIR$ IVDEP
            DO k = data%ISTADA( ig ), data%ISTADA( ig1 ) - 1
              ll = data%ICNA( k )

!  include contributions from the first n variables only

              IF ( ll <= n ) THEN
                IF ( JTRANS ) THEN
                  CJAC( ll, icon ) = gi * data%A( k )
                ELSE
                  CJAC( icon, ll ) = gi * data%A( k )
                END IF
              END IF
            END DO
          END IF
        END DO
      END IF
      status = 0
      RETURN

!  unsuccessful returns

  930 CONTINUE
      IF ( data%out > 0 ) WRITE( data%out,                                     &
        "( ' ** SUBROUTINE CCFG: error flag raised during SIF evaluation' )" )
      status = 3
      RETURN

!  non-executable statements

 2000 FORMAT( ' ** SUBROUTINE CCFG: Increase the leading dimension of CJAC' ) 
 2010 FORMAT( ' ** SUBROUTINE CCFG: Increase the second dimension of CJAC' )

!  end of subroutine CCFG

      END SUBROUTINE CCFG

