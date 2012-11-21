!-*-*-*-*-*-  CUTEst_interfaces m O D U l E  *-*-*-*-*-*-*-*

!  Nick Gould, for CUTEst productions
!  Copyright reserved
!  October 21st 2012

   Module CUTEst_interfaces

!  Interface blocks for CUTE fortran tools
        
!!$    USE CUTEst_precis
       USE CUTEST
     Implicit None

     Interface

!  Interface block for unconstrained tools

       Subroutine UDIMEN( input, n )
         Use CUTEst_precis
         Integer, Intent( IN ) :: input
         Integer, Intent( OUT ) :: n
       End Subroutine UDIMEN

       Subroutine USETUP( data, input, iout, n, X, BL, BU, nmax )
         Use CUTEst_precis
         TYPE ( CUTEST_data_type ) :: data
         Integer, Intent( IN ) :: input, iout, nmax
         Integer, Intent( OUT ) :: n
         Real ( KIND = wp ), Intent( OUT ), Dimension( nmax ) :: X, BL, BU
       End Subroutine USETUP

       Subroutine UNAMES( data, n, PNAME, VNAME )
         Use CUTEst_precis
         TYPE ( CUTEST_data_type ) :: data
         Integer, Intent( IN ) :: n
         Character * 10, Intent( OUT ) :: PNAME
         Character * 10, Intent( OUT ), Dimension( n ) :: VNAME
       End Subroutine UNAMES

       Subroutine PBNAME( data, n, PNAME )
         Use CUTEst_precis
         TYPE ( CUTEST_data_type ) :: data
         Integer, Intent( IN ) :: n
         Character * 10, Intent( OUT ) :: PNAME
       End Subroutine PBNAME

       Subroutine VARNAMES( data, n, VNAME )
         Use CUTEst_precis
         TYPE ( CUTEST_data_type ) :: data
         Integer, Intent( IN ) :: n
         Character * 10, Intent( OUT ), Dimension( n ) :: VNAME
       End Subroutine VARNAMES

       Subroutine UVARTY( data, n, IVARTY )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) :: n
         Integer, Intent( OUT ) :: IVARTY( n )
       End Subroutine UVARTY

       Subroutine UFN( data, n, X, f )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) ::  n
         Real ( KIND = wp ), Intent( OUT ) :: f
         Real ( KIND = wp ), Intent( IN ), Dimension( n ) :: X
       End Subroutine UFN

       Subroutine UGR( data, n, X, G )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) ::  n
         Real ( KIND = wp ), Intent( IN ), Dimension( n ) :: X
         Real ( KIND = wp ), Intent( OUT ), Dimension( n ) :: G
       End Subroutine UGR

       Subroutine UOFG( data, n, X, f, G, GRAD )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) ::  n
         Real ( KIND = wp ), Intent( OUT ) :: f
         Real ( KIND = wp ), Intent( IN ), Dimension( n ) :: X
         Real ( KIND = wp ), Intent( OUT ), Dimension( n ) :: G
         Logical, Intent( IN ) :: GRAD
       End Subroutine UOFG

       Subroutine UDH( data, n, X, lh1, H )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) ::  n, lh1
         Real ( KIND = wp ), Intent( IN ), Dimension( n ) :: X
         Real ( KIND = wp ), Intent( OUT ), Dimension( lh1, n ) :: H
       End Subroutine UDH

       Subroutine UGRDH( data, n, X, G, lh1, H )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) :: n, lh1
         Real ( KIND = wp ), Intent( IN ), Dimension( n ) :: X
         Real ( KIND = wp ), Intent( OUT ), Dimension( n ) :: G
         Real ( KIND = wp ), Intent( OUT ), Dimension( lh1, n ) :: H
       End Subroutine UGRDH

       Subroutine UDIMSH( data, nnzh )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( OUT ) :: nnzh
       End Subroutine UDIMSH

       Subroutine USH( data, n, X, nnzh, lh, H, IRNH, ICNH )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) :: n, lh
         Integer, Intent( OUT ) :: nnzh
         Integer, Intent( OUT ), Dimension( lh ) :: IRNH, ICNH
         Real ( KIND = wp ), Intent( IN ), Dimension( n ) :: X
         Real ( KIND = wp ), Intent( OUT ), Dimension( lh ) :: H
       End Subroutine USH

       Subroutine UDIMSE( data, ne, nnzh, nzirnh )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( OUT ) :: ne, nnzh, nzirnh
       End Subroutine UDIMSE

       Subroutine UEH( data, n, X, ne, IRNHI, lirnhi, le,                     &
                       IPRNHI, HI, lhi, IPRHI, BYROWS )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) :: n, le, lirnhi, lhi 
         Integer, Intent( IN ) :: ne
         Logical, Intent( IN ) :: BYROWS
         Integer, Intent( OUT ), Dimension( lirnhi ) :: IRNHI
         Integer, Intent( OUT ), Dimension( le ) :: IPRNHI, IPRHI
         Real ( KIND = wp ), Intent( IN ), Dimension( n ) :: X
         Real ( KIND = wp ), Intent( OUT ), Dimension( lhi ) :: HI
       End Subroutine UEH

       Subroutine UGRSH( data, n, X, G, nnzh, lh, H, IRNH, ICNH )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) :: n, lh
         Integer, Intent( OUT ) :: nnzh
         Integer, Intent( OUT ), Dimension( lh ) :: IRNH, ICNH
         Real ( KIND = wp ), Intent( IN ), Dimension( n ) :: X
         Real ( KIND = wp ), Intent( OUT ), Dimension( n ) :: G
         Real ( KIND = wp ), Intent( OUT ), Dimension( lh ) :: H
       End Subroutine UGRSH

       Subroutine UGREH( data, n, X, G, ne, IRNHI, lirnhi, le,                &
                         IPRNHI, HI, lhi, IPRHI, BYROWS )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) :: n, le, lirnhi, lhi 
         Integer, Intent( OUT ) :: ne
         Logical, Intent( IN ) :: BYROWS
         Integer, Intent( OUT ), Dimension( lirnhi ) :: IRNHI
         Integer, Intent( OUT ), Dimension( le ) :: IPRNHI, IPRHI
         Real ( KIND = wp ), Intent( IN ), Dimension( n ) :: X
         Real ( KIND = wp ), Intent( OUT ), Dimension( n ) :: G
         Real ( KIND = wp ), Intent( OUT ), Dimension( lhi ) :: HI
       End Subroutine UGREH

       Subroutine UPROD( data, n, GOTH, X, P, Result )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) :: n
         Logical, Intent( IN ) :: GOTH
         Real ( KIND = wp ), Intent( IN ), Dimension( n ) :: X, P
         Real ( KIND = wp ), Intent( OUT ), Dimension( n ) :: Result
       End Subroutine UPROD

       Subroutine UBANDH( data, n, GOTH, X, nsemib, BANDH, lbandh )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) ::  n, nsemib, lbandh
         Logical, Intent( IN ) ::  GOTH
         Real ( KIND = wp ), Intent( IN ), Dimension( n ) ::  X
         Real ( KIND = wp ), Intent( OUT ),                                    &
              Dimension( 0 : lbandh, n ) ::  BANDH
       End Subroutine UBANDH

       !  Interface block for constrained tools

       Subroutine CDIMEN( input, n, m )
         Use CUTEst_precis
         Integer, Intent( IN ) :: input
         Integer, Intent( OUT ) :: n, m
       End Subroutine CDIMEN

       Subroutine CSETUP( data, input, iout, n, m, X, BL, BU, nmax, EQUATN,    &
                          LINEAR, V, CL, CU, mmax, EFIRST, LFIRST, NVFRST )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) ::  input, iout, nmax, mmax
         Integer, Intent( OUT ) ::  n, m
         Logical, Intent( IN ) ::  EFIRST, LFIRST, NVFRST
         Real ( KIND = wp ), Intent( OUT ), Dimension( nmax ) :: X, BL, BU
         Real ( KIND = wp ), Intent( OUT ), Dimension( mmax ) :: V, CL, CU
         Logical, Intent( OUT ), Dimension( mmax ) :: EQUATN, LINEAR
       End Subroutine CSETUP

       Subroutine CNAMES( data, n, m, PNAME, VNAME, GNAME )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) :: n, m
         Character * 10, Intent( OUT ) :: PNAME
         Character * 10, Intent( OUT ), Dimension( n ) :: VNAME
         Character * 10, Intent( OUT ), Dimension( m ) :: GNAME
       End Subroutine CNAMES

       Subroutine CONNAMES( data, m, GNAME )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) :: m
         Character * 10, Intent( OUT ), Dimension( m ) :: GNAME
       End Subroutine CONNAMES

       Subroutine CVARTY( data, n, IVARTY )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) :: n
         Integer, Intent( OUT ) :: IVARTY( n )
       End Subroutine CVARTY

       Subroutine CFN( data, n, m, X, f, lc, C )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) :: n, m, lc
         Real ( KIND = wp ), Intent( OUT ) :: f
         Real ( KIND = wp ), Intent( IN ), Dimension( n ) :: X
         Real ( KIND = wp ), Intent( OUT ), Dimension( lc ) :: C
       End Subroutine CFN

       Subroutine CGR( data, n, m, X, GRLAGF, lv, V, G, JTRANS,               &
                       lcjac1, lcjac2, CJAC )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) :: n, m, lv, lcjac1, lcjac2
         Logical, Intent( IN ) :: GRLAGF, JTRANS
         Real ( KIND = wp ), Intent( IN ), Dimension( n ) :: X
         Real ( KIND = wp ), Intent( IN ), Dimension( lv ) :: V
         Real ( KIND = wp ), Intent( OUT ), Dimension( n ) :: G
         Real ( KIND = wp ), Intent( OUT ),                                    &
              Dimension( lcjac1, lcjac2 ) :: CJAC
       End Subroutine CGR

       Subroutine COFG( data, n, X, f, G, GRAD )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) :: n
         Real ( KIND = wp ), Intent( OUT ) :: f
         Real ( KIND = wp ), Intent( IN ), Dimension( n ) :: X
         Real ( KIND = wp ), Intent( OUT ), Dimension( n ) :: G
         Logical, Intent( IN ) :: GRAD
       End Subroutine COFG

       Subroutine CDIMSJ( data, nnzj )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( OUT ) :: nnzj
       End Subroutine CDIMSJ

       Subroutine CSGR( data, n, m, GRLAGF, lv, V, X, nnzj,                    &
                        lcjac, CJAC, INDVAR, INDFUN )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) :: n, m, lv, lcjac
         Integer, Intent( OUT ) :: nnzj
         Logical, Intent( IN ) :: GRLAGF
         Integer, Intent( OUT ), Dimension( lcjac ) :: INDVAR, INDFUN
         Real ( KIND = wp ), Intent( IN ), Dimension( n ) :: X
         Real ( KIND = wp ), Intent( IN ), Dimension( lv ) :: V
         Real ( KIND = wp ), Intent( OUT ), Dimension( lcjac ) :: CJAC
       End Subroutine CSGR

       Subroutine CCFG( data, n, m, X, lc, C, JTRANS, lcjac1, lcjac2, CJAC,   &
                        GRAD )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) :: n, m, lc, lcjac1, lcjac2
         Real ( KIND = wp ), Intent( IN ), Dimension( n ) :: X
         Real ( KIND = wp ), Intent( OUT ), Dimension( lc ) :: C
         Real ( KIND = wp ), Intent( OUT ),                                   &
              Dimension( lcjac1, lcjac2 ) :: CJAC
         Logical, Intent( IN ) :: JTRANS, GRAD
       End Subroutine CCFG

       Subroutine CSCFG( data, n, m, X, lc, C, nnzj, lcjac, CJAC,             &
                         INDVAR, INDFUN, GRAD )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) :: n, m, lc, lcjac
         Integer, Intent( OUT ) :: nnzj
         Logical, Intent( IN ) :: GRAD
         Integer, Intent( OUT ), Dimension( lcjac ) :: INDVAR, INDFUN
         Real ( KIND = wp ), Intent( IN ), Dimension( n ) :: X
         Real ( KIND = wp ), Intent( OUT ), Dimension( lc ) :: C
         Real ( KIND = wp ), Intent( OUT ), Dimension( lcjac ) :: CJAC
       End Subroutine CSCFG

       Subroutine CCIFG( data, n, icon, X, ci, GCI, GRAD )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) ::  n, icon
         Logical, Intent( IN ) :: GRAD
         Real ( KIND = wp ), Intent( OUT ) :: ci
         Real ( KIND = wp ), Intent( IN ), Dimension( n ) :: X
         Real ( KIND = wp ), Intent( OUT ), Dimension( n ) :: GCI
       End Subroutine CCIFG

       Subroutine CSCIFG( data, n, icon, X, ci, nnzgci, lgci, GCI, INDVAR,    &
                          GRAD )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) :: n, icon, lgci
         Integer, Intent( OUT ) :: nnzgci
         Logical, Intent( IN ) :: GRAD
         Integer, Intent( OUT ), Dimension( lgci ) :: INDVAR
         Real ( KIND = wp ), Intent( OUT ) :: ci
         Real ( KIND = wp ), Intent( IN ), Dimension( n ) :: X
         Real ( KIND = wp ), Intent( OUT ), Dimension( lgci ) :: GCI
       End Subroutine CSCIFG

       Subroutine CDH( data, n, m, X, lv, V, lh1, H )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) ::  n, m, lv, lh1
         Real ( KIND = wp ), Intent( IN ), Dimension( n ) :: X
         Real ( KIND = wp ), Intent( IN ), Dimension( lv ) :: V
         Real ( KIND = wp ), Intent( OUT ), Dimension( lh1, n ) :: H
       End Subroutine CDH

       Subroutine CIDH( data, n, X, iprob, lh1, H )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) ::  n, iprob, lh1
         Real ( KIND = wp ), Intent( IN ), Dimension( n ) :: X
         Real ( KIND = wp ), Intent( OUT ), Dimension( lh1, n ) :: H
       End Subroutine CIDH

       Subroutine CGRDH( data, n, m, X, GRLAGF, lv, V, G,                     &
                         JTRANS, lcjac1, lcjac2, CJAC, lh1, H )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) :: n, m, lv, lh1, lcjac1, lcjac2
         Logical, Intent( IN ) :: GRLAGF, JTRANS
         Real ( KIND = wp ), Intent( IN ), Dimension( n ) :: X
         Real ( KIND = wp ), Intent( OUT ), Dimension( n ) :: G
         Real ( KIND = wp ), Intent( IN ), Dimension( lv ) :: V
         Real ( KIND = wp ), Intent( OUT ), Dimension( lh1, n ) :: H
         Real ( KIND = wp ), Intent( OUT ),                                   &
              Dimension( lcjac1, lcjac2 ) :: CJAC
       End Subroutine CGRDH

       Subroutine CDIMSH( data, nnzh )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( OUT ) :: nnzh
       End Subroutine CDIMSH

       Subroutine CSH( data, n, m, X, lv, V, nnzh, lh, H, IRNH, ICNH )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) :: n, m, lv, lh
         Integer, Intent( OUT ) :: nnzh
         Integer, Intent( OUT ), Dimension( lh ) :: IRNH, ICNH
         Real ( KIND = wp ), Intent( IN ), Dimension( n ) :: X
         Real ( KIND = wp ), Intent( IN ), Dimension( lv ) :: V
         Real ( KIND = wp ), Intent( OUT ), Dimension( lh ) :: H
       End Subroutine CSH

       Subroutine CSH1( data, n, m, X, lv, V, nnzh, lh, H, IRNH, ICNH )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) :: n, m, lv, lh
         Integer, Intent( OUT ) :: nnzh
         Integer, Intent( OUT ), Dimension( lh ) :: IRNH, ICNH
         Real ( KIND = wp ), Intent( IN ), Dimension( n ) :: X
         Real ( KIND = wp ), Intent( IN ), Dimension( lv ) :: V
         Real ( KIND = wp ), Intent( OUT ), Dimension( lh ) :: H
       End Subroutine CSH1

       Subroutine CISH( data, n, X, iprob, nnzh, lh, H, IRNH, ICNH )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) :: n, iprob, lh
         Integer, Intent( OUT ) :: nnzh
         Integer, Intent( OUT ), Dimension( lh ) :: IRNH, ICNH
         Real ( KIND = wp ), Intent( IN ), Dimension( n ) :: X
         Real ( KIND = wp ), Intent( OUT ), Dimension( lh ) :: H
       End Subroutine CISH

       Subroutine CDIMSE( data, ne, nnzh, nzirnh )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( OUT ) :: ne, nnzh, nzirnh
       End Subroutine CDIMSE

       Subroutine CEH( data, n, m, X, lv, V, ne, IRNHI, lirnhi, le,           &
                       IPRNHI, HI, lhi, IPRHI, BYROWS )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) :: n, m, lv, le, lirnhi, lhi 
         Integer, Intent( OUT ) :: ne
         Logical, Intent( IN ) :: BYROWS
         Integer, Intent( OUT ), Dimension( lirnhi ) :: IRNHI
         Integer, Intent( OUT ), Dimension( le ) :: IPRNHI, IPRHI
         Real ( KIND = wp ), Intent( IN ), Dimension( n ) :: X
         Real ( KIND = wp ), Intent( IN ), Dimension( lv ) :: V
         Real ( KIND = wp ), Intent( OUT ), Dimension( lhi ) :: HI
       End Subroutine CEH

       Subroutine CSGRSH( data, n, m, X, GRLAGF, lv, V, nnzj, lcjac, CJAC,    &
                          INDVAR, INDFUN, nnzh, lh, H, IRNH, ICNH )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) ::  n, m, lv, lcjac, lh
         Integer, Intent( OUT ) :: nnzj, nnzh
         Logical, Intent( IN ) ::  GRLAGF
         Integer, Intent( OUT ), Dimension( lcjac ) :: INDVAR, INDFUN
         Integer, Intent( OUT ), Dimension( lh ) :: IRNH, ICNH
         Real ( KIND = wp ), Intent( IN ), Dimension( n ) :: X
         Real ( KIND = wp ), Intent( IN ), Dimension( lv ) :: V
         Real ( KIND = wp ), Intent( OUT ), Dimension( lh ) :: H
         Real ( KIND = wp ), Intent( OUT ), Dimension( lcjac ) :: CJAC
       End Subroutine CSGRSH

       Subroutine CSGREH( data, n, m, X, GRLAGF, lv, V, nnzj, lcjac, CJAC,    &
                          INDVAR, INDFUN, ne, IRNHI, lirnhi, le, IPRNHI,      &
                          HI, lhi, IPRHI, BYROWS )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) :: n, m, lv, lcjac, le, lirnhi, lhi 
         Integer, Intent( OUT ) :: ne, nnzj
         Logical, Intent( IN ) :: GRLAGF, BYROWS
         Integer, Intent( IN ), Dimension( lcjac ) :: INDVAR, INDFUN
         Integer, Intent( OUT ), Dimension( lirnhi ) :: IRNHI
         Integer, Intent( OUT ), Dimension( le ) :: IPRNHI, IPRHI
         Real ( KIND = wp ), Intent( IN ), Dimension( n ) :: X
         Real ( KIND = wp ), Intent( IN ), Dimension( lv ) :: V
         Real ( KIND = wp ), Intent( OUT ), Dimension( lhi ) :: HI
         Real ( KIND = wp ), Intent( OUT ), Dimension( lcjac ) :: CJAC
       End Subroutine CSGREH

       Subroutine CPROD( data, n, m, GOTH, X, lv, V, P, Result )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) :: n, m, lv
         Logical, Intent( IN ) :: GOTH
         Real ( KIND = wp ), Intent( IN ), Dimension( n ) :: X, P
         Real ( KIND = wp ), Intent( IN ), Dimension( lv ) :: V
         Real ( KIND = wp ), Intent( OUT ), Dimension( n ) :: Result
       End Subroutine CPROD
   
       Subroutine CPROD1( data, n, m, GOTH, X, lv, V, P, Result )
         TYPE ( CUTEST_data_type ) :: data
         Use CUTEst_precis
         Integer, Intent( IN ) :: n, m, lv
         Logical, Intent( IN ) :: GOTH
         Real ( KIND = wp ), Intent( IN ), Dimension( n ) :: X, P
         Real ( KIND = wp ), Intent( IN ), Dimension( lv ) :: V
         Real ( KIND = wp ), Intent( OUT ), Dimension( n ) :: Result
       End Subroutine CPROD1
   
     End Interface
         
!  End of module CUTEst_interfaces

   End Module CUTEst_interfaces
      
