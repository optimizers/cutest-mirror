! ( Last modified on Mon Feb 25 15:03:37 CST 2002 )
      SUBROUTINE CSCFG ( data, n, m, X, lc, C, nnzj, lcjac, CJAC, &
                         INDVAR, INDFUN, GRAD )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n, m, lc, nnzj, lcjac
      LOGICAL :: GRAD
      INTEGER :: INDVAR( lcjac ), INDFUN( lcjac )
      REAL ( KIND = wp ) :: X( n ), C( lc ), CJAC ( lcjac )

!  *********************************************************************
!  *                                                                   *
!  *             OBSOLETE TOOL! Replaced by CCFSG                     *
!  *                                                                   *
!  *********************************************************************


      WRITE( iout, 1000 )
      CALL CCFSG ( n, m, X, lc, C, nnzj, lcjac, CJAC, &
                   INDVAR, INDFUN, GRAD )
      RETURN

!  Non executable statement

 1000 FORMAT( ' ** SUBROUTINE CSCFG: this tool is obsolete! ', &
              ' Please use CCFSG instead.' )

!  end of CSCFG.

      END

