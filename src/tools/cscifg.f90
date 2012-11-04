! ( Last modified on Mon Feb 25 15:03:37 CST 2002 )
      SUBROUTINE CSCIFG ( data, n, icon, X, ci, nnzgci, lgci, GCI, &
                          INDVAR, GRAD )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: n, icon, nnzgci, lgci
      LOGICAL :: GRAD
      INTEGER :: INDVAR( lgci )
      REAL ( KIND = wp ) :: ci
      REAL ( KIND = wp ) :: X( n ), GCI( lgci )

!  *********************************************************************
!  *                                                                   *
!  *             OBSOLETE TOOL! Replaced by CCIFSG                     *
!  *                                                                   *
!  *********************************************************************


      WRITE( iout, 1000 )
      CALL CCIFSG( n, icon, X, ci, nnzgci, lgci, GCI, &
                   INDVAR, GRAD )
      RETURN

!  Non executable statement

 1000 FORMAT( ' ** SUBROUTINE CSCIFG: this tool is obsolete! ', &
              ' Please use CCIFSG instead.' )

!  end of CSCIFG.

      END

