! ( Last modified on Mon Feb 25 15:03:37 CST 2002 )
      SUBROUTINE CSCIFG ( data, N, ICON, X, CI, NNZGCI, LGCI, GCI, &
                          INDVAR, GRAD )
      USE CUTEST
      TYPE ( CUTEST_data_type ) :: data
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER :: N, ICON, NNZGCI, LGCI
      LOGICAL :: GRAD
      INTEGER :: INDVAR( LGCI )
      REAL ( KIND = wp ) :: CI
      REAL ( KIND = wp ) :: X( N ), GCI( LGCI )

!  *********************************************************************
!  *                                                                   *
!  *             OBSOLETE TOOL! Replaced by CCIFSG                     *
!  *                                                                   *
!  *********************************************************************


      WRITE( IOUT, 1000 )
      CALL CCIFSG( N, ICON, X, CI, NNZGCI, LGCI, GCI, &
                   INDVAR, GRAD )
      RETURN

!  Non executable statement

 1000 FORMAT( ' ** SUBROUTINE CSCIFG: this tool is obsolete! ', &
              ' Please use CCIFSG instead.' )

!  end of CSCIFG.

      END

