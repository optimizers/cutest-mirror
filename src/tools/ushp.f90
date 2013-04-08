! THIS VERSION: CUTEST 1.0 - 07/04/2013 AT 10:00 GMT.

!-*-*-*-*-*-*-*-  C U T E S T    U S H P  S U B R O U T I N E  -*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 2003 version released in CUTEst, 7th April 2013

      SUBROUTINE CUTEST_ushp( status, n, nnzh, lh, H_row, H_col )
      USE CUTEST

!  dummy arguments

      INTEGER, INTENT( IN ) :: n, lh
      INTEGER, INTENT( OUT ) :: nnzh, status
      INTEGER, INTENT( OUT ), DIMENSION( lh ) :: H_row, H_col

!  ---------------------------------------------------------------
!  compute the spasity pattern of the Hessian matrix of a group 
!  partially separable function. The upper triangle of the Hessian 
!  is stored in coordinate form, i.e., the entry has row index 
!  H_row(i) and column index H_col(i) for i = 1, ...., nnzh
!  ---------------------------------------------------------------

      CALL CUTEST_ushp_threadsafe( CUTEST_data_global,                         &
                                   CUTEST_work_global( 1 ), status, n,         &
                                   nnzh, lh, H_row, H_col )
      RETURN

!  end of subroutine CUTEST_ushp

      END SUBROUTINE CUTEST_ushp

!-*-*-  C U T E S T   U S H P _ t h r e a d s a f e   S U B R O U T I N E  -*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 2003 version released in CUTEst, 7th April 2013

      SUBROUTINE CUTEST_ushp_threadsafe( data, work, status, n,                &
                                         nnzh, lh, H_row, H_col )
      USE CUTEST

!  dummy arguments

      TYPE ( CUTEST_data_type ), INTENT( IN ) :: data
      TYPE ( CUTEST_work_type ), INTENT( INOUT ) :: work
      INTEGER, INTENT( IN ) :: n, lh
      INTEGER, INTENT( OUT ) :: nnzh, status
      INTEGER, INTENT( OUT ), DIMENSION( lh ) :: H_row, H_col

!  ---------------------------------------------------------------
!  compute the spasity pattern of the Hessian matrix of a group 
!  partially separable function. The upper triangle of the Hessian 
!  is stored in coordinate form, i.e., the entry has row index 
!  H_row(i) and column index H_col(i) for i = 1, ...., nnzh
!  ---------------------------------------------------------------

!  local variables

      INTEGER :: alloc_status
      CHARACTER ( LEN = 80 ) :: bad_alloc = REPEAT( ' ', 80 )

!  determine the Hessian pattern

      CALL CUTEST_assemble_hessian_pattern(                                    &
             n, data%ng, data%nel, data%ntotel, data%nvrels, data%nvargp,      &
             data%ISTADH, data%IELVAR, data%IELING, data%ISTADG,               &
             data%ISTAEV, data%ISTAGV, data%ISVGRP, data%GXEQX,                &
             0, data%out, data%out, work%io_buffer,                            &
             n, status, alloc_status, bad_alloc,                               &
             work%array_status, work%lh_row, work%lh_col,                      &
             work%H_row, work%H_col,                                           &
             work%LINK_col, work%POS_in_H, work%llink, work%lpos, nnzh )

!  check for errors

      IF ( status > 0 ) RETURN

!  record the sparse Hessian

      H_row( : nnzh ) = work%H_row( : nnzh )
      H_col( : nnzh ) = work%H_col( : nnzh )

!  update the counters for the report tool

      work%nc2oh = work%nc2oh + 1
      status = 0
      RETURN

!  end of subroutine CUTEST_ushp_threadsafe

      END SUBROUTINE CUTEST_ushp_threadsafe
