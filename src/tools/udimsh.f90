! THIS VERSION: CUTEST 1.0 - 28/11/2012 AT 13:30 GMT.

!-*-*-*-*-*-*-  C U T E S T   U D I M S H    S U B R O U T I N E  -*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 2003 version released in CUTEst, 28th December 2012

      SUBROUTINE CUTEST_udimsh( status, nnzh )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      INTEGER, INTENT( OUT ) :: status, nnzh

!  -----------------------------------------------------------------------
!  compute the space required to store the Hessian matrix of the objective 
!  function of a problem initially written in Standard Input Format (SIF)

!  NB. USETUP must have been called first

!  The upper triangle of the Hessian is stored in coordinate form,
!  i.e., the entry H_val(i) has row index H_row(i) for i = 1, ...., nnzh
!  -----------------------------------------------------------------------

      CALL CUTEST_udimsh_threadsafe( CUTEST_data_global,                       &
                                     CUTEST_work_global( 1 ), status, nnzh )
      RETURN

!  end of subroutine CUTEST_udimsh

      END SUBROUTINE CUTEST_udimsh

!-*-  C U T E S T   U D I M S H _ t h r e a d s a f e   S U B R O U T I N E  -*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 77 version originally released in CUTEr, August 1999
!   fortran 2003 version released in CUTEst, 24th November 2012

      SUBROUTINE CUTEST_udimsh_threadsafe( data, work, status, nnzh )
      USE CUTEST
      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  dummy arguments

      TYPE ( CUTEST_data_type ), INTENT( IN ) :: data
      TYPE ( CUTEST_work_type ), INTENT( INOUT ) :: work
      INTEGER, INTENT( OUT ) :: status, nnzh

!  -----------------------------------------------------------------------
!  compute the space required to store the Hessian matrix of the objective 
!  function of a problem initially written in Standard Input Format (SIF)

!  NB. USETUP must have been called first

!  The upper triangle of the Hessian is stored in coordinate form,
!  i.e., the entry H_val(i) has row index H_row(i) for i = 1, ...., nnzh
!  -----------------------------------------------------------------------

!  local variables

      INTEGER :: lirnh, alloc_status
      CHARACTER ( LEN = 80 ) :: bad_alloc = REPEAT( ' ', 80 )

!  guess required space

      lirnh = 10000
      CALL CUTEST_size_sparse_hessian(                                         &
                      data%n, data%ng, data%nel, data%ntotel, data%nvrels,     &
                      data%nvargp, data%IELVAR, data%IELING,                   &
                      data%ISTADG, data%ISTAEV, data%ISTAGV, data%ISVGRP,      &
                      data%GXEQX, data%out, work%io_buffer, status,            &
                      alloc_status, bad_alloc, lirnh, work%H_row,              &
                      work%LINK_col, work%POS_in_H, work%llink, work%lpos,     &
                      nnzh )

      RETURN

!  end of subroutine CUTEST_udimsh_threadsafe

      END SUBROUTINE CUTEST_udimsh_threadsafe
