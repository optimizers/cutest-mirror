# Digital/Compaq f95 compiler
#
#  Fortran compilation and loading
#

FORTRAN='f95'
BASIC='-c -check nopower'
LIBCMD=''
MODCMD='-module $MOD -I$MOD'
MVMODS=':'
OPTIMIZATION='-O -arch host -tune host'
NOOPTIMIZATION='-O0'
DEBUG=
OPENMP=
F77='-fixed'
F90='-w'
F95='-w'
NOFMAIN='-nofor_main'
CCONDEF=
USUAL=
SPECIAL=
F77SUFFIX=f
F95SUFFIX=f90
TIMER=GEN
BLAS=
LAPACK=
HSL=
METIS='-lgalahad_metis'
PARDISO='-lgalahad_pardiso'
WSMP='-lgalahad_wsmp'
NOT95=IS95
NOT64=IS64
BINSHELL=sh
