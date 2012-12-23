! THIS VERSION: CUTEST 1.0 - 22/11/2012 AT 08:45 GMT.

!-*-*-*-*-*-*-  C U T E S T    C S C I F G    S U B R O U T I N E  -*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!   fortran 77 version originally released in CUTE, December 1990
!   replaced in CUTEr, December 1990, by subroutine CCIFSG
!   fortran 2003 version released in CUTEst, 22nd November 2012

      SUBROUTINE CUTEST_cscifg( data, work, status, n, icon, X, ci,                  &
                                nnzgci, lgci, GCI_val, GCI_var, grad )
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

!  *********************************************************************
!  *                                                                   *
!  *             OBSOLETE TOOL! Replaced by CCIFSG                     *
!  *                                                                   *
!  *********************************************************************

      IF ( data%out > 0 ) WRITE( data%out, "( ' ** SUBROUTINE CSCIFG: this ',  &
    &  'tool is obsolete!  Please use CCIFSG instead.' )" )

      CALL CUTEST_ccifsg( data, work, status, n, icon, X, ci,                        &
                          nnzgci, lgci, GCI_val, GCI_var, grad )
      RETURN

!  end of subroutine CUTEST_cscifg

      END SUBROUTINE CUTEST_cscifg

