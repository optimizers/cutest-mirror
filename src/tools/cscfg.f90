! THIS VERSION: CUTEST 1.0 - 22/11/2012 AT 08:45 GMT.

!-*-*-*-*-*-*-  C U T E S T    C S C F G    S U B R O U T I N E  -*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!   fortran 77 version originally released in CUTE, December 1990
!   replaced in CUTEr, December 1990, by subroutine CCFSG
!   fortran 2003 version released in CUTEst, 22nd November 2012

      SUBROUTINE CSCFG( data, status, n, m, X, C, nnzj, lcjac, CJAC,           &
                        INDVAR, INDFUN, grad )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      TYPE ( CUTEST_data_type ), INTENT( INOUT ) :: data
      INTEGER, INTENT( IN ) :: n, m, lcjac
      INTEGER, INTENT( OUT ) :: status, nnzj
      LOGICAL, INTENT( IN ) :: grad
      INTEGER, INTENT( OUT ), DIMENSION( lcjac ) :: INDVAR, INDFUN
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: X
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( m ) :: C
      REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( lcjac ) :: CJAC

!  *********************************************************************
!  *                                                                   *
!  *             OBSOLETE TOOL! Replaced by CCFSG                      *
!  *                                                                   *
!  *********************************************************************

      IF ( data%out > 0 ) WRITE( data%out, "( ' ** SUBROUTINE CSCIFG: this ',  &
    &  'tool is obsolete! Please use CCFSG instead.' )" )
      CALL CCFSG( data, status, n, m, X, C, nnzj, lcjac, CJAC, INDVAR, INDFUN, &
                  grad )
      RETURN

!  end of subroutine CSCFG

      END SUBROUTINE CSCFG

