/*
 * ======================================================================
 *
 * cutest.h
 * Data type definitions, constants definitions and function prototypes
 * to interface the CUTEst testing environment Fortran library with C
 *
 * This header file is built from different sources and different authors
 * have contributed to it. In any case, many thanks to Andreas Waechter.
 *
 * D. Orban. CUTEr version, July 10 2002.
 * Nick Gould, CUTEst evolution, January 4 2013.
 *             Boolean logicals provided, August 21 2013
 * ======================================================================
 */

#include <stdlib.h>
#include <stdbool.h>

#ifndef CUTEST_DOT_H_INCLUDED
#define CUTEST_DOT_H_INCLUDED
#endif

/*
 * Define name of main() function on a
 * compiler by compiler basis.
 */
#ifdef Isg95
#define MAINENTRY MAIN_
#endif
#ifdef Ispgf
#define MAINENTRY MAIN_
#endif
#ifdef Isifr
#define MAINENTRY MAIN__
#endif

#ifndef MAINENTRY
#define MAINENTRY main
#endif

/*
 * Define Fortran types for integer and double precision
 * The following choices are from f2c.h
 */
/* typedef long int integer; */
typedef int      integer;
typedef float    real;
typedef double   doublereal;
typedef _Bool    logical;
/* typedef bool    logical; */
#define FALSE_ (0)     /* Fortran FALSE */
#define TRUE_  (1)     /* Fortran  TRUE */
/* #define max(a,b) ((a)>(b)?(a):(b)) */

#define ZERO     0e0
#define ONE      1e0
#define CUTE_INF 1e20    /* 'infinity' in CUTEst interface */
#define FSTRING_LEN 10   /* Length of Fortran strings     */

/* AIX does not append underscore to Fortran subroutine names */
#ifdef _AIX
#define FUNDERSCORE(a)   a
#else
#define FUNDERSCORE(a)   a##_
#endif

typedef struct VarTypes {
	int nbnds, neq, nlin, nrange, nlower, nupper,
		nineq, nineq_lin, nineq_nlin, neq_lin,
		neq_nlin;
} VarTypes;

/*
 * Define shortcuts for the CUTEst library functions,
 * and try to avoid the trailing underscore.
 *
 */

#define CUTEST_usetup   FUNDERSCORE(cutest_usetup)
#define CUTEST_csetup   FUNDERSCORE(cutest_cint_csetup)

#define CUTEST_udimen   FUNDERSCORE(cutest_udimen)
#define CUTEST_udimsh   FUNDERSCORE(cutest_udimsh)
#define CUTEST_udimse   FUNDERSCORE(cutest_udimse)
#define CUTEST_uvartype FUNDERSCORE(cutest_uvartype)
#define CUTEST_unames   FUNDERSCORE(cutest_unames)
#define CUTEST_ureport  FUNDERSCORE(cutest_ureport)

#define CUTEST_cdimen   FUNDERSCORE(cutest_cdimen)
#define CUTEST_cdimsj   FUNDERSCORE(cutest_cdimsj)
#define CUTEST_cdimsh   FUNDERSCORE(cutest_cdimsh)
#define CUTEST_cdimse   FUNDERSCORE(cutest_cdimse)
#define CUTEST_cdstats  FUNDERSCORE(cutest_cstats)
#define CUTEST_cvartype FUNDERSCORE(cutest_cvartype)
#define CUTEST_cnames   FUNDERSCORE(cutest_cnames)
#define CUTEST_creport  FUNDERSCORE(cutest_creport)

#define CUTEST_connames FUNDERSCORE(cutest_connames)
#define CUTEST_pname    FUNDERSCORE(cutest_pname)
#define CUTEST_probname FUNDERSCORE(cutest_probname)
#define CUTEST_varnames FUNDERSCORE(cutest_varnames)

#define CUTEST_ufn      FUNDERSCORE(cutest_ufn)
#define CUTEST_ugr      FUNDERSCORE(cutest_ugr)
#define CUTEST_uofg     FUNDERSCORE(cutest_cint_uofg)
#define CUTEST_ubandh   FUNDERSCORE(cutest_ubandh)
#define CUTEST_udh      FUNDERSCORE(cutest_udh)
#define CUTEST_ush      FUNDERSCORE(cutest_ush)
#define CUTEST_ueh      FUNDERSCORE(cutest_cint_ueh)
#define CUTEST_ugrdh    FUNDERSCORE(cutest_ugrdh)
#define CUTEST_ugrsh    FUNDERSCORE(cutest_ugrsh)
#define CUTEST_ugreh    FUNDERSCORE(cutest_cint_ugreh)
#define CUTEST_uhprod   FUNDERSCORE(cutest_cint_uhprod)

#define CUTEST_cfn      FUNDERSCORE(cutest_cfn)
#define CUTEST_cofg     FUNDERSCORE(cutest_cint_cofg)
#define CUTEST_cofsg    FUNDERSCORE(cutest_cint_cofsg)
#define CUTEST_ccfg     FUNDERSCORE(cutest_cint_ccfg)
#define CUTEST_clfg     FUNDERSCORE(cutest_cint_clfg)
#define CUTEST_cgr      FUNDERSCORE(cutest_cint_cgr)
#define CUTEST_csgr     FUNDERSCORE(cutest_cint_csgr)
#define CUTEST_ccfsg    FUNDERSCORE(cutest_cint_ccfsg)
#define CUTEST_ccifg    FUNDERSCORE(cutest_cint_ccifg)
#define CUTEST_ccifsg   FUNDERSCORE(cutest_cint_ccifsg)
#define CUTEST_cgrdh    FUNDERSCORE(cutest_cint_cgrdh)
#define CUTEST_cdh      FUNDERSCORE(cutest_cdh)
#define CUTEST_csh      FUNDERSCORE(cutest_csh)
#define CUTEST_cshc     FUNDERSCORE(cutest_cshc)
#define CUTEST_ceh      FUNDERSCORE(cutest_cint_ceh)
#define CUTEST_cidh     FUNDERSCORE(cutest_cidh)
#define CUTEST_cish     FUNDERSCORE(cutest_cish)
#define CUTEST_csgrsh   FUNDERSCORE(cutest_cint_csgrsh)
#define CUTEST_csgreh   FUNDERSCORE(cutest_cint_csgreh)
#define CUTEST_chprod   FUNDERSCORE(cutest_cint_chprod)
#define CUTEST_chcprod  FUNDERSCORE(cutest_cint_chcprod)
#define CUTEST_cjprod   FUNDERSCORE(cutest_cint_cjprod)

#define CUTEST_uterminate FUNDERSCORE(cutest_uterminate)
#define CUTEST_cterminate FUNDERSCORE(cutest_cterminate)

#define FORTRAN_open  FUNDERSCORE(fortran_open)
#define FORTRAN_close FUNDERSCORE(fortran_close)

/*
 * Prototypes for CUTEst FORTRAN routines found in libcutest.a
 * See  http://ccpforge.cse.rl.ac.uk/gf/project/cutest/
 */

/* Setup routines */
void CUTEST_usetup( integer *status, integer *funit, integer *iout, 
              integer *io_buffer, integer *n, doublereal *x,
	      doublereal *bl, doublereal *bu );
void CUTEST_csetup( integer *status, integer *funit, integer *iout, 
             integer *io_buffer, integer *n, integer *m,
	      doublereal *x, doublereal *bl, doublereal *bu, 
              doublereal *v, doublereal *cl, doublereal *cu, 
	      logical *equatn, logical *linear, 
              integer *e_order, integer *l_order, integer *v_order );

/* Unconstrained dimensioning and report routines */
void CUTEST_udimen( integer *status, integer *funit, integer *n );
void CUTEST_udimsh( integer *status, integer *nnzh );
void CUTEST_udimse( integer *status, integer *ne, integer *nzh,
                    integer *nzirnh );
void CUTEST_uvartype( integer *status, integer *n, integer *ivarty );
void CUTEST_unames( integer *status, integer *n, char *pname, char *vnames );
void CUTEST_ureport( integer *status, doublereal *calls, doublereal *time );

/* Constrained dimensioning and report routines */
void CUTEST_cdimen( integer *status, integer *funit, integer *n, integer *m );
void CUTEST_cdimsj( integer *status, integer *nnzj );
void CUTEST_cdimsh( integer *status, integer *nnzh );
void CUTEST_cdimse( integer *status, integer *ne, integer *nzh, 
                    integer *nzirnh );
void CUTEST_cstats( integer *status, integer *nonlinear_variables_objective, 
                    integer *nonlinear_variables_constraints,
                    integer *equality_constraints, 
                    integer *linear_constraints );
void CUTEST_cvartype( integer *status, integer *n, integer *ivarty );
void CUTEST_cnames( integer *status, integer *n, integer *m, char *pname, 
                    char *vnames, char *gnames );
void CUTEST_creport( integer *status, doublereal *calls, doublereal *time );

void CUTEST_connames( integer *status, integer *m, char *gname );
void CUTEST_pname( integer *status, integer *funit, char *pname );
void CUTEST_probname( integer *status, char *pname );
void CUTEST_varnames( integer *status, integer *n, char *vname );

/* Unconstrained optimization routines */
void CUTEST_ufn( integer *status, integer *n, doublereal *x, doublereal *f );
void CUTEST_ugr( integer *status, integer *n, doublereal *x, doublereal *g );
void CUTEST_uofg( integer *status, integer *n, doublereal *x, doublereal *f, 
                  doublereal *g, logical *grad );
void CUTEST_ubandh( integer *status, integer *n, doublereal *x, 
                    integer *nsemib, doublereal *bandh, integer *lbandh, 
                    integer *maxsbw );
void CUTEST_udh( integer *status, integer *n, doublereal *x, integer *lh1, 
          doublereal *h );
void CUTEST_ush( integer *status, integer *n, doublereal *x, integer *nnzh, 
           integer *lh, doublereal *h, integer *irnh, integer *icnh );
void CUTEST_ueh( integer *status, integer *n, doublereal *x, integer *ne, 
          integer *le, integer *iprnhi, integer *iprhi, integer *lirnhi, 
          integer *irnhi, integer *lhi, doublereal *hi, logical *byrows );

void CUTEST_ugrdh( integer *status, integer *n, doublereal *x, doublereal *g, 
            integer *lh1, doublereal *h);
void CUTEST_ugrsh( integer *status, integer *n, doublereal *x, doublereal *g, 
             integer *nnzh,
	     integer *lh, doublereal *h, integer *irnh, integer *icnh );
void CUTEST_ugreh( integer *status, integer *n, doublereal *x, doublereal *g, 
             integer *ne,
             integer *le, integer *iprnhi, integer *iprhi, integer *lirnhi, 
             integer *irnhi, integer *lhi, doublereal *hi, logical *byrows );
void CUTEST_uhprod( integer *status, integer *n, logical *goth, doublereal *x, 
             doublereal *p, doublereal *q );

/* Constrained optimization routines */
void CUTEST_cfn( integer *status,  integer *n, integer *m, doublereal *x, 
          doublereal *f, doublereal *c );
void CUTEST_cofg( integer *status, integer *n, doublereal *x, doublereal *f, 
           doublereal *g, logical *grad );
void CUTEST_cofsg( integer *status, integer *n, doublereal *x, doublereal *f, 
                   integer *nnzg, integer *lg, 
                   doublereal *sg, integer *ivsg, logical *grad );
void CUTEST_ccfg( integer *status, integer *n, integer *m, doublereal *x, 
	    doublereal *c, logical *jtrans, integer *lcjac1, integer *lcjac2,
	    doublereal *cjac, logical *grad );
void CUTEST_clfg( integer *status, integer *n, integer *m, doublereal *x, 
           doublereal *v, doublereal *f, doublereal *g, logical *grad );
void CUTEST_cgr( integer *status,  integer *n, integer *m, doublereal *x, 
            doublereal *v, logical *grlagf, doublereal *g, logical *jtrans,
	    integer *lcjac1, integer *lcjac2, doublereal *cjac );
void CUTEST_csgr( integer *status, integer *n, integer *m, 
            doublereal *x, doublereal *v, logical *grlagf, 
            integer *nnzj, integer *lcjac,
	    doublereal *cjac, integer *indvar, integer *indfun );
void CUTEST_ccfsg( integer *status,  integer *n, integer *m, doublereal *x, 
              doublereal *c, integer *nnzj, integer *lcjac,
	      doublereal *cjac, integer *indvar, integer *indfun,
	      logical *grad );
void CUTEST_ccifg( integer *status,  integer *n, integer *i, doublereal *x, 
              doublereal *ci, doublereal *gci, logical *grad );
void CUTEST_ccifsg( integer *status, integer *n, integer *i, doublereal *x, 
              doublereal *ci, integer *nnzsgc, integer *lsgci, 
              doublereal *sgci, integer *ivsgci, logical *grad );
void CUTEST_cgrdh( integer *status, integer *n, integer *m, doublereal *x, 
             doublereal *v, logical *grlagf, doublereal *g, logical *jtrans,
	     integer *lcjac1, integer *lcjac2, doublereal *cjac,
	     integer *lh1, doublereal *h );
void CUTEST_cdh( integer *status, integer *n, integer *m, doublereal *x, 
           doublereal *v, integer *lh1, doublereal *h );
void CUTEST_csh( integer *status, integer *n, integer *m, doublereal *x, 
           doublereal *v, integer *nnzh, integer *lh, doublereal *h, 
           integer *irnh, integer *icnh );
void CUTEST_cshc( integer *status, integer *n, integer *m, doublereal *x, 
           doublereal *v, integer *nnzh, integer *lh, doublereal *h, 
           integer *irnh, integer *icnh );
void CUTEST_ceh( integer *status, integer *n, integer *m, doublereal *x, 
           integer *lv, doublereal *v, integer *ne, 
           integer *le, integer *iprnhi, integer *iprhi, integer *lirnhi, 
           integer *irnhi, integer *lhi, doublereal *hi, logical *byrows );
void CUTEST_cidh( integer *status, integer *n, doublereal *x, integer *iprob, 
            integer *lh1, doublereal *h );
void CUTEST_cish( integer *status, integer *n, doublereal *x, integer *iprob, 
            integer *nnzh,
	    integer *lh, doublereal *h, integer *irnh, integer *icnh );
void CUTEST_csgrsh( integer *status, integer *n, integer *m, doublereal *x, 
	      doublereal *v, logical *grlagf, integer *nnzj, integer *lcjac,
	      doublereal *cjac, integer *indvar, integer *indfun,
	      integer *nnzh, integer *lh, doublereal *h, integer *irnh,
	      integer *icnh );
void CUTEST_csgreh( integer *status, integer *n, integer *m, doublereal *x, 
	      doublereal *v, logical *grlagf, integer *nnzj, integer *lcjac,
	      doublereal *cjac, integer *indvar, integer *indfun,
	      integer *ne, 
              integer *le, integer *iprnhi, integer *iprhi, integer *lirnhi, 
              integer *irnhi, integer *lhi, doublereal *hi, logical *byrows );
void CUTEST_chprod( integer *status, integer *n, integer *m, logical *goth, 
              doublereal *x, doublereal *v, doublereal *p, doublereal *q );
void CUTEST_chcprod( integer *status, integer *n, integer *m, logical *goth, 
               doublereal *x, doublereal *v, doublereal *p, doublereal *q );
void CUTEST_cjprod( integer *status, integer *n, integer *m, logical *gotj, 
              logical *jtrans, doublereal *x, doublereal *p, integer *lp, 
              doublereal *r, integer *lr );

/* Termination routines */
void CUTEST_uterminate( integer *status );
void CUTEST_cterminate( integer *status );

/* FORTRAN auxiliary subroutines to retrieve stream unit numbers */
void FORTRAN_open(  integer *funit, char *fname, integer *ierr );
void FORTRAN_close( integer *funit, integer *ierr );

/*
 * Memory allocation shortcuts
 */

void *CUTEst_malloc( void *object, int length, size_t s );
void *CUTEst_calloc( void *object, int length, size_t s );
void *CUTEst_realloc( void *object, int length, size_t s );
void  CUTEst_free( void **object );

#ifndef MALLOC
#define MALLOC(object,length,type)  object = (type *)CUTEst_malloc(object,length,sizeof(type))
#endif
#ifndef CALLOC
#define CALLOC(object,length,type)  object = (type *)CUTEst_calloc(object,length,sizeof(type))
#endif
#ifndef REALLOC
#define REALLOC(object,length,type) object = (type *)CUTEst_realloc(object,length,sizeof(type))
#endif
#ifndef FREE
#define FREE(object) CUTEst_free((void **)(&(object)))
#endif

