! THIS VERSION: CUTEST 1.0 - 26/11/2012 AT 14:20 GMT.

!-*-*-*-*-*-*-  C U T E S T   C D I M S E    S U B R O U T I N E  -*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   fortran 77 version originally released in CUTEr, November 1994
!   fortran 2003 version released in CUTEst, 26th November 2012

      SUBROUTINE CDIMSE( data, status, ne, he_val_ne, he_row_ne )
      USE CUTEST

!  dummy arguments

      TYPE ( CUTEST_data_type ), INTENT( IN ) :: data
      INTEGER, INTENT( OUT ) :: status, ne, he_val_ne, he_row_ne

!  --------------------------------------------------------------------------
!  compute the number of elements and the space required to store the Hessian 
!  matrix of the Lagrangian function of a problem initially written in 
!  Standard Input Format (SIF)

!  The matrix is represented in "finite element format", i.e., 

!           ne
!      H = sum H_e, 
!          e=1

!  where each element H_i involves a small subset of the rows of H.
!  H is stored as a list of the row indices involved in each element
!  and the upper triangle of H_e (stored by rows or columns). 

!  ne (integer) number of elements
!  he_val_ne (integer) number of entries needed to store the real values of H. 
!         Specifically, the sum of the number of entries in the upper triangle
!         of each H_e
!  he_row_ne (integer) number of entries needed to store the integer entries of
!         H. Specifically, the sum of the row dimensions of each H_e
!  ---------------------------------------------------------------------------

!  History -
!   fortran 77 version released in CUTEr as CDIMSE, November 25th 1994
!   fortran 2003 version released in CUTEst, 26th November 2012

      CALL CUTEST_size_element_hessian( data%ng, data%ISTADG, data%ISTAGV,     &
                                        data%GXEQX, ne, he_val_ne, he_row_ne,  &
                                        status )
      RETURN

!  end of subroutine CDIMSE

      END SUBROUTINE CDIMSE
