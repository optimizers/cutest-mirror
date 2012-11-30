! THIS VERSION: CUTEST 1.0 - 25/11/2012 AT 13:35 GMT.

!-*-*-*-*-*-*-  C U T E S T   C D I M S H    S U B R O U T I N E  -*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 77 version originally released in CUTEr, August 1999
!   fortran 2003 version released in CUTEst, 24th November 2012

      SUBROUTINE CDIMSH( data, status, nnzh )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      TYPE ( CUTEST_data_type ), INTENT( INOUT ) :: data
      INTEGER, INTENT( OUT ) :: status, nnzh

!  ------------------------------------------------------------------------
!  Compute the space required to store the Hessian matrix of the Lagrangian 
!  function of a problem initially written in Standard Input Format (SIF)

!  NB. CSETUP must have been called first

!  The upper triangle of the Hessian is stored in coordinate form,
!  i.e., the entry H_val(i) has row index H_row(i) for i = 1, ...., nnzh.

!  ------------------------------------------------------------------------

!  local variables

      INTEGER :: lirnh, alloc_status
      CHARACTER ( LEN = 80 ) :: bad_alloc = REPEAT( ' ', 80 )

!  guess required space

      lirnh = 10000
      CALL CUTEST_size_sparse_hessian(                                         &
                      data%n, data%ng, data%nel, data%ntotel, data%nvrels,     &
                      data%nvargp, data%IELVAR, data%IELING,                   &
                      data%ISTADG, data%ISTAEV, data%ISTAGV, data%ISVGRP,      &
                      data%GXEQX, data%out, data%io_buffer, status,            &
                      alloc_status, bad_alloc, lirnh, data%H_row,              &
                      data%LINK_col, data%POS_in_H, data%llink, data%lpos,     &
                      nnzh )

      RETURN

!  end of subroutine CDIMSH

      END SUBROUTINE CDIMSH
