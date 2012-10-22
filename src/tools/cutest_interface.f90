!-*-*-*-*-*-  CUTEst_interfaces M O D U L E  *-*-*-*-*-*-*-*

!  Nick Gould, for CUTEst productions
!  Copyright reserved
!  October 21st 2012

   Module CUTEst_interfaces

!  Interface blocks for CUTE fortran tools
        
!!$      USE CUTEst_precis
      Implicit None

      Interface

!  Interface block for unconstrained tools

         Subroutine UDIMEN( INPUT, N )
           Use CUTEst_precis
           Integer, Intent( IN ) :: INPUT
           Integer, Intent( OUT ) :: N
         End Subroutine UDIMEN

         Subroutine USETUP( INPUT , IOUT  , N , X , BL, BU, NMAX )
           Use CUTEst_precis
           Integer, Intent( IN ) :: INPUT, IOUT, NMAX
           Integer, Intent( OUT ) :: N
           Real ( KIND = wp ), Intent( OUT ), Dimension( NMAX ) :: X, BL, BU
         End Subroutine USETUP

         Subroutine UNAMES( N, PNAME, VNAME )
           Use CUTEst_precis
           Integer, Intent( IN ) :: N
           Character * 10, Intent( OUT ) :: PNAME
           Character * 10, Intent( OUT ), Dimension( N ) :: VNAME
         End Subroutine UNAMES

         Subroutine PBNAME( N, PNAME )
           Use CUTEst_precis
           Integer, Intent( IN ) :: N
           Character * 10, Intent( OUT ) :: PNAME
         End Subroutine PBNAME

         Subroutine VARNAMES( N, VNAME )
           Use CUTEst_precis
           Integer, Intent( IN ) :: N
           Character * 10, Intent( OUT ), Dimension( N ) :: VNAME
         End Subroutine VARNAMES

         Subroutine UVARTY( N, IVARTY )
           Use CUTEst_precis
           Integer, Intent( IN ) :: N
           Integer, Intent( OUT ) :: IVARTY( N )
         End Subroutine UVARTY

         Subroutine UFN( N, X, F )
           Use CUTEst_precis
           Integer, Intent( IN ) ::  N
           Real ( KIND = wp ), Intent( OUT ) :: F
           Real ( KIND = wp ), Intent( IN ), Dimension( N ) :: X
         End Subroutine UFN

         Subroutine UGR( N, X, G )
           Use CUTEst_precis
           Integer, Intent( IN ) ::  N
           Real ( KIND = wp ), Intent( IN ), Dimension( N ) :: X
           Real ( KIND = wp ), Intent( OUT ), Dimension( N ) :: G
         End Subroutine UGR

         Subroutine UOFG( N, X, F, G, GRAD )
           Use CUTEst_precis
           Integer, Intent( IN ) ::  N
           Real ( KIND = wp ), Intent( OUT ) :: F
           Real ( KIND = wp ), Intent( IN ), Dimension( N ) :: X
           Real ( KIND = wp ), Intent( OUT ), Dimension( N ) :: G
           Logical, Intent( IN ) :: GRAD
         End Subroutine UOFG

         Subroutine UDH( N, X, LH1, H )
           Use CUTEst_precis
           Integer, Intent( IN ) ::  N, LH1
           Real ( KIND = wp ), Intent( IN ), Dimension( N ) :: X
           Real ( KIND = wp ), Intent( OUT ), Dimension( LH1, N ) :: H
         End Subroutine UDH

         Subroutine UGRDH( N, X, G, LH1, H )
           Use CUTEst_precis
           Integer, Intent( IN ) :: N, LH1
           Real ( KIND = wp ), Intent( IN ), Dimension( N ) :: X
           Real ( KIND = wp ), Intent( OUT ), Dimension( N ) :: G
           Real ( KIND = wp ), Intent( OUT ), Dimension( LH1, N ) :: H
         End Subroutine UGRDH

         Subroutine UDIMSH( NNZH )
           Use CUTEst_precis
           Integer, Intent( OUT ) :: NNZH
         End Subroutine UDIMSH

         Subroutine USH( N, X, NNZH, LH, H, IRNH, ICNH )
           Use CUTEst_precis
           Integer, Intent( IN ) :: N, LH
           Integer, Intent( OUT ) :: NNZH
           Integer, Intent( OUT ), Dimension( LH ) :: IRNH, ICNH
           Real ( KIND = wp ), Intent( IN ), Dimension( N ) :: X
           Real ( KIND = wp ), Intent( OUT ), Dimension( LH ) :: H
         End Subroutine USH

         Subroutine UDIMSE( NE, NNZH, NZIRNH )
           Use CUTEst_precis
           Integer, Intent( OUT ) :: NE, NNZH, NZIRNH
         End Subroutine UDIMSE

         Subroutine UEH( N , X , NE    , IRNHI , LIRNHI, LE    ,              &
              IPRNHI, HI    , LHI   , IPRHI , BYROWS )
           Use CUTEst_precis
           Integer, Intent( IN ) :: N , LE    , LIRNHI, LHI 
           Integer, Intent( IN ) :: NE
           Logical, Intent( IN ) :: BYROWS
           Integer, Intent( OUT ), Dimension( LIRNHI ) :: IRNHI
           Integer, Intent( OUT ), Dimension( LE ) :: IPRNHI, IPRHI
           Real ( KIND = wp ), Intent( IN ), Dimension( N ) :: X
           Real ( KIND = wp ), Intent( OUT ), Dimension( LHI ) :: HI
         End Subroutine UEH

         Subroutine UGRSH( N, X, G, NNZH, LH, H, IRNH, ICNH )
           Use CUTEst_precis
           Integer, Intent( IN ) :: N, LH
           Integer, Intent( OUT ) :: NNZH
           Integer, Intent( OUT ), Dimension( LH ) :: IRNH, ICNH
           Real ( KIND = wp ), Intent( IN ), Dimension( N ) :: X
           Real ( KIND = wp ), Intent( OUT ), Dimension( N ) :: G
           Real ( KIND = wp ), Intent( OUT ), Dimension( LH ) :: H
         End Subroutine UGRSH

         Subroutine UGREH( N , X , G, NE , IRNHI , LIRNHI, LE    ,            &
              IPRNHI, HI    , LHI   , IPRHI , BYROWS )
           Use CUTEst_precis
           Integer, Intent( IN ) :: N, LE, LIRNHI, LHI 
           Integer, Intent( OUT ) :: NE
           Logical, Intent( IN ) :: BYROWS
           Integer, Intent( OUT ), Dimension( LIRNHI ) :: IRNHI
           Integer, Intent( OUT ), Dimension( LE ) :: IPRNHI, IPRHI
           Real ( KIND = wp ), Intent( IN ), Dimension( N ) :: X
           Real ( KIND = wp ), Intent( OUT ), Dimension( N ) :: G
           Real ( KIND = wp ), Intent( OUT ), Dimension( LHI ) :: HI
         End Subroutine UGREH

         Subroutine UPROD( N, GOTH, X, P, Result )
           Use CUTEst_precis
           Integer, Intent( IN ) :: N
           Logical, Intent( IN ) :: GOTH
           Real ( KIND = wp ), Intent( IN ), Dimension( N ) :: X, P
           Real ( KIND = wp ), Intent( OUT ), Dimension( N ) :: Result
         End Subroutine UPROD

         Subroutine UBANDH( N, GOTH, X, NSEMIB, BANDH, LBANDH )
           Use CUTEst_precis
           Integer, Intent( IN ) ::  N, NSEMIB, LBANDH
           Logical, Intent( IN ) ::  GOTH
           Real ( KIND = wp ), Intent( IN ), Dimension( N ) ::  X
           Real ( KIND = wp ), Intent( OUT ),                              &
                Dimension( 0 : LBANDH ) ::  BANDH
         End Subroutine UBANDH

         !  Interface block for constrained tools

         Subroutine CDIMEN( INPUT, N, M )
           Use CUTEst_precis
           Integer, Intent( IN ) :: INPUT
           Integer, Intent( OUT ) :: N, M
         End Subroutine CDIMEN

         Subroutine CSETUP( INPUT , IOUT  , N , M , X , BL , BU   , NMAX,     &
              EQUATN, LINEAR, V , CL, CU     , MMAX , EFIRST,   &
              LFIRST, NVFRST )
           Use CUTEst_precis
           Integer, Intent( IN ) ::  INPUT , IOUT  , NMAX  , MMAX
           Integer, Intent( OUT ) ::  N, M
           Logical, Intent( IN ) ::  EFIRST, LFIRST, NVFRST
           Real ( KIND = wp ), Intent( OUT ), Dimension( NMAX ) :: X, BL, BU
           Real ( KIND = wp ), Intent( OUT ), Dimension( MMAX ) :: V, CL, CU
           Logical, Intent( OUT ), Dimension( MMAX ) :: EQUATN, LINEAR
         End Subroutine CSETUP

         Subroutine CNAMES( N, M, PNAME, VNAME, GNAME )
           Use CUTEst_precis
           Integer, Intent( IN ) :: N, M
           Character * 10, Intent( OUT ) :: PNAME
           Character * 10, Intent( OUT ), Dimension( N ) :: VNAME
           Character * 10, Intent( OUT ), Dimension( M ) :: GNAME
         End Subroutine CNAMES

         Subroutine CONNAMES( M, GNAME )
           Use CUTEst_precis
           Integer, Intent( IN ) :: M
           Character * 10, Intent( OUT ), Dimension( M ) :: GNAME
         End Subroutine CONNAMES

         Subroutine CVARTY( N, IVARTY )
           Use CUTEst_precis
           Integer, Intent( IN ) :: N
           Integer, Intent( OUT ) :: IVARTY( N )
         End Subroutine CVARTY

         Subroutine CFN( N , M , X , F , LC, C )
           Use CUTEst_precis
           Integer, Intent( IN ) :: N , M , LC
           Real ( KIND = wp ), Intent( OUT ) :: F
           Real ( KIND = wp ), Intent( IN ), Dimension( N ) :: X
           Real ( KIND = wp ), Intent( OUT ), Dimension( LC ) :: C
         End Subroutine CFN

         Subroutine CGR( N , M , X     , GRLAGF, LV, V , G     , JTRANS,      &
              LCJAC1, LCJAC2, CJAC  )
           Use CUTEst_precis
           Integer, Intent( IN ) :: N , M , LV    , LCJAC1, LCJAC2
           Logical, Intent( IN ) :: GRLAGF, JTRANS
           Real ( KIND = wp ), Intent( IN ), Dimension( N ) :: X
           Real ( KIND = wp ), Intent( IN ), Dimension( LV ) :: V
           Real ( KIND = wp ), Intent( OUT ), Dimension( N ) :: G
           Real ( KIND = wp ), Intent( OUT ),                              &
                Dimension( LCJAC1, LCJAC2 ) :: CJAC
         End Subroutine CGR

         Subroutine COFG( N, X, F, G, GRAD )
           Use CUTEst_precis
           Integer, Intent( IN ) :: N
           Real ( KIND = wp ), Intent( OUT ) :: F
           Real ( KIND = wp ), Intent( IN ), Dimension( N ) :: X
           Real ( KIND = wp ), Intent( OUT ), Dimension( N ) :: G
           Logical, Intent( IN ) :: GRAD
         End Subroutine COFG

         Subroutine CDIMSJ( NNZJ )
           Use CUTEst_precis
           Integer, Intent( OUT ) :: NNZJ
         End Subroutine CDIMSJ

         Subroutine CSGR( N , M , GRLAGF, LV, V , X     , NNZJ  ,             &
              LCJAC , CJAC  , INDVAR, INDFUN )
           Use CUTEst_precis
           Integer, Intent( IN ) :: N , M , LV, LCJAC
           Integer, Intent( OUT ) :: NNZJ
           Logical, Intent( IN ) :: GRLAGF
           Integer, Intent( OUT ), Dimension( LCJAC ) :: INDVAR, INDFUN
           Real ( KIND = wp ), Intent( IN ), Dimension( N ) :: X
           Real ( KIND = wp ), Intent( IN ), Dimension( LV ) :: V
           Real ( KIND = wp ), Intent( OUT ), Dimension( LCJAC ) :: CJAC
         End Subroutine CSGR

         Subroutine CCFG( N, M, X, LC, C, JTRANS, LCJAC1, LCJAC2, CJAC, GRAD )
           Use CUTEst_precis
           Integer, Intent( IN ) :: N, M, LC, LCJAC1, LCJAC2
           Real ( KIND = wp ), Intent( IN ), Dimension( N ) :: X
           Real ( KIND = wp ), Intent( OUT ), Dimension( LC ) :: C
           Real ( KIND = wp ), Intent( OUT ),                              &
                Dimension( LCJAC1, LCJAC2 ) :: CJAC
           Logical, Intent( IN ) :: JTRANS, GRAD
         End Subroutine CCFG

         Subroutine CSCFG( N, M, X, LC, C, NNZJ, LCJAC, CJAC,                 &
              INDVAR, INDFUN, GRAD )
           Use CUTEst_precis
           Integer, Intent( IN ) :: N, M, LC, LCJAC
           Integer, Intent( OUT ) :: NNZJ
           Logical, Intent( IN ) :: GRAD
           Integer, Intent( OUT ), Dimension( LCJAC ) :: INDVAR, INDFUN
           Real ( KIND = wp ), Intent( IN ), Dimension( N ) :: X
           Real ( KIND = wp ), Intent( OUT ), Dimension( LC ) :: C
           Real ( KIND = wp ), Intent( OUT ), Dimension( LCJAC ) :: CJAC
         End Subroutine CSCFG

         Subroutine CCIFG( N, ICON, X, CI, GCI, GRAD )
           Use CUTEst_precis
           Integer, Intent( IN ) ::  N, ICON
           Logical, Intent( IN ) :: GRAD
           Real ( KIND = wp ), Intent( OUT ) :: CI
           Real ( KIND = wp ), Intent( IN ), Dimension( N ) :: X
           Real ( KIND = wp ), Intent( OUT ), Dimension( N ) :: GCI
         End Subroutine CCIFG

         Subroutine CSCIFG( N, ICON, X, CI, NNZGCI, LGCI, GCI, INDVAR, GRAD )
           Use CUTEst_precis
           Integer, Intent( IN ) :: N, ICON, LGCI
           Integer, Intent( OUT ) :: NNZGCI
           Logical, Intent( IN ) :: GRAD
           Integer, Intent( OUT ), Dimension( LGCI ) :: INDVAR
           Real ( KIND = wp ), Intent( OUT ) :: CI
           Real ( KIND = wp ), Intent( IN ), Dimension( N ) :: X
           Real ( KIND = wp ), Intent( OUT ), Dimension( LGCI ) :: GCI
         End Subroutine CSCIFG

         Subroutine CDH( N, M, X, LV, V, LH1, H )
           Use CUTEst_precis
           Integer, Intent( IN ) ::  N, M, LV, LH1
           Real ( KIND = wp ), Intent( IN ), Dimension( N ) :: X
           Real ( KIND = wp ), Intent( IN ), Dimension( LV ) :: V
           Real ( KIND = wp ), Intent( OUT ), Dimension( LH1, N ) :: H
         End Subroutine CDH

         Subroutine CIDH( N, X, IPROB, LH1, H )
           Use CUTEst_precis
           Integer, Intent( IN ) ::  N, IPROB, LH1
           Real ( KIND = wp ), Intent( IN ), Dimension( N ) :: X
           Real ( KIND = wp ), Intent( OUT ), Dimension( LH1, N ) :: H
         End Subroutine CIDH

         Subroutine CGRDH( N , M , X     , GRLAGF, LV , V, G     ,            &
              JTRANS, LCJAC1, LCJAC2, CJAC  , LH1, H     )
           Use CUTEst_precis
           Integer, Intent( IN ) :: N , M , LV    , LH1   , LCJAC1, LCJAC2
           Logical, Intent( IN ) :: GRLAGF, JTRANS
           Real ( KIND = wp ), Intent( IN ), Dimension( N ) :: X
           Real ( KIND = wp ), Intent( OUT ), Dimension( N ) :: G
           Real ( KIND = wp ), Intent( IN ), Dimension( LV ) :: V
           Real ( KIND = wp ), Intent( OUT ), Dimension( LH1, N ) :: H
           Real ( KIND = wp ), Intent( OUT ),                              &
                Dimension( LCJAC1, LCJAC2 ) :: CJAC
         End Subroutine CGRDH

         Subroutine CDIMSH( NNZH )
           Use CUTEst_precis
           Integer, Intent( OUT ) :: NNZH
         End Subroutine CDIMSH

         Subroutine CSH( N, M, X, LV, V, NNZH, LH, H, IRNH, ICNH  )
           Use CUTEst_precis
           Integer, Intent( IN ) :: N, M, LV, LH
           Integer, Intent( OUT ) :: NNZH
           Integer, Intent( OUT ), Dimension( LH ) :: IRNH, ICNH
           Real ( KIND = wp ), Intent( IN ), Dimension( N ) :: X
           Real ( KIND = wp ), Intent( IN ), Dimension( LV ) :: V
           Real ( KIND = wp ), Intent( OUT ), Dimension( LH ) :: H
         End Subroutine CSH

         Subroutine CSH1( N, M, X, LV, V, NNZH, LH, H, IRNH, ICNH  )
           Use CUTEst_precis
           Integer, Intent( IN ) :: N, M, LV, LH
           Integer, Intent( OUT ) :: NNZH
           Integer, Intent( OUT ), Dimension( LH ) :: IRNH, ICNH
           Real ( KIND = wp ), Intent( IN ), Dimension( N ) :: X
           Real ( KIND = wp ), Intent( IN ), Dimension( LV ) :: V
           Real ( KIND = wp ), Intent( OUT ), Dimension( LH ) :: H
         End Subroutine CSH1

         Subroutine CISH( N, X, IPROB, NNZH, LH, H, IRNH, ICNH  )
           Use CUTEst_precis
           Integer, Intent( IN ) :: N, IPROB, LH
           Integer, Intent( OUT ) :: NNZH
           Integer, Intent( OUT ), Dimension( LH ) :: IRNH, ICNH
           Real ( KIND = wp ), Intent( IN ), Dimension( N ) :: X
           Real ( KIND = wp ), Intent( OUT ), Dimension( LH ) :: H
         End Subroutine CISH

         Subroutine CDIMSE( NE, NNZH, NZIRNH )
           Use CUTEst_precis
           Integer, Intent( OUT ) :: NE, NNZH, NZIRNH
         End Subroutine CDIMSE

         Subroutine CEH( N , M , X , LV, V , NE, IRNHI , LIRNHI, LE    ,      &
              IPRNHI, HI    , LHI   , IPRHI , BYROWS )
           Use CUTEst_precis
           Integer, Intent( IN ) :: N, M, LV, LE, LIRNHI, LHI 
           Integer, Intent( OUT ) :: NE
           Logical, Intent( IN ) :: BYROWS
           Integer, Intent( OUT ), Dimension( LIRNHI ) :: IRNHI
           Integer, Intent( OUT ), Dimension( LE ) :: IPRNHI, IPRHI
           Real ( KIND = wp ), Intent( IN ), Dimension( N ) :: X
           Real ( KIND = wp ), Intent( IN ), Dimension( LV ) :: V
           Real ( KIND = wp ), Intent( OUT ), Dimension( LHI ) :: HI
         End Subroutine CEH

         Subroutine CSGRSH( N , M , X     , GRLAGF, LV, V , NNZJ  ,           &
              LCJAC , CJAC  , INDVAR, INDFUN, NNZH  ,           &
              LH, H , IRNH  , ICNH  )
           Use CUTEst_precis
           Integer, Intent( IN ) ::  N, M, LV, LCJAC , LH
           Integer, Intent( OUT ) :: NNZJ, NNZH
           Logical, Intent( IN ) ::  GRLAGF
           Integer, Intent( OUT ), Dimension( LCJAC ) :: INDVAR, INDFUN
           Integer, Intent( OUT ), Dimension( LH ) :: IRNH, ICNH
           Real ( KIND = wp ), Intent( IN ), Dimension( N ) :: X
           Real ( KIND = wp ), Intent( IN ), Dimension( LV ) :: V
           Real ( KIND = wp ), Intent( OUT ), Dimension( LH ) :: H
           Real ( KIND = wp ), Intent( OUT ), Dimension( LCJAC ) :: CJAC
         End Subroutine CSGRSH

         Subroutine CSGREH( N , M , X     , GRLAGF, LV, V , NNZJ  , LCJAC ,   &
              CJAC  , INDVAR, INDFUN, NE    , IRNHI , LIRNHI,   &
              LE    , IPRNHI, HI    , LHI   , IPRHI , BYROWS )
           Use CUTEst_precis
           Integer, Intent( IN ) :: N, M, LV, LCJAC, LE, LIRNHI, LHI 
           Integer, Intent( OUT ) :: NE, NNZJ
           Logical, Intent( IN ) :: GRLAGF, BYROWS
           Integer, Intent( IN ), Dimension( LCJAC ) :: INDVAR, INDFUN
           Integer, Intent( OUT ), Dimension( LIRNHI ) :: IRNHI
           Integer, Intent( OUT ), Dimension( LE ) :: IPRNHI, IPRHI
           Real ( KIND = wp ), Intent( IN ), Dimension( N ) :: X
           Real ( KIND = wp ), Intent( IN ), Dimension( LV ) :: V
           Real ( KIND = wp ), Intent( OUT ), Dimension( LHI ) :: HI
           Real ( KIND = wp ), Intent( OUT ), Dimension( LCJAC ) :: CJAC
         End Subroutine CSGREH

         Subroutine CPROD( N , M , GOTH  , X , LV, V , P , Result )
           Use CUTEst_precis
           Integer, Intent( IN ) :: N , M , LV
           Logical, Intent( IN ) :: GOTH
           Real ( KIND = wp ), Intent( IN ), Dimension( N ) :: X, P
           Real ( KIND = wp ), Intent( IN ), Dimension( LV ) :: V
           Real ( KIND = wp ), Intent( OUT ), Dimension( N ) :: Result
         End Subroutine CPROD
   
         Subroutine CPROD1( N , M , GOTH  , X , LV, V , P , Result )
           Use CUTEst_precis
           Integer, Intent( IN ) :: N , M , LV
           Logical, Intent( IN ) :: GOTH
           Real ( KIND = wp ), Intent( IN ), Dimension( N ) :: X, P
           Real ( KIND = wp ), Intent( IN ), Dimension( LV ) :: V
           Real ( KIND = wp ), Intent( OUT ), Dimension( N ) :: Result
         End Subroutine CPROD1
   
      End Interface
         
!  End of module CUTEst_interfaces

   End Module CUTEst_interfaces
      
