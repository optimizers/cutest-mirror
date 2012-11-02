! ( Last modified on Mon Feb 25 15:03:37 CST 2002 )
      SUBROUTINE CSCFG ( data, N, M, X, LC, C, NNZJ, LCJAC, CJAC, &
                         INDVAR, INDFUN, GRAD )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: N, M, LC, NNZJ, LCJAC
      LOGICAL :: GRAD
      INTEGER :: INDVAR( LCJAC ), INDFUN( LCJAC )
      REAL ( KIND = wp ) :: X( N ), C( LC ), CJAC ( LCJAC )

!  *********************************************************************
!  *                                                                   *
!  *             OBSOLETE TOOL! Replaced by CCFSG                     *
!  *                                                                   *
!  *********************************************************************


      WRITE( IOUT, 1000 )
      CALL CCFSG ( N, M, X, LC, C, NNZJ, LCJAC, CJAC, &
                   INDVAR, INDFUN, GRAD )
      RETURN

!  Non executable statement

 1000 FORMAT( ' ** SUBROUTINE CSCFG: this tool is obsolete! ', &
              ' Please use CCFSG instead.' )

!  end of CSCFG.

      END

