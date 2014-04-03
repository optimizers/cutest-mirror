/* ====================================================
 * CUTEst interface for cg_descent     April. 2, 2014
 *
 * W. Hager
 *
 * (Based on CUTEr gencma.c of D. Orban, Feb 3, 2003)
 * (CUTEst evolution, Nick Gould, Apr 2, 2014)
 * ====================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define CG_DESCENTMA

#ifdef __cplusplus
extern "C" {   /* To prevent C++ compilers from mangling symbols */
#endif

/*
#include "../../include/cuter.h"
#include "../../include/cg_user.h"
*/
#include "cutest.h"
#include "cg_user.h"

#ifdef Isg95
#define MAINENTRY MAIN_
#else
#define MAINENTRY main
#endif

/* prototypes */
double cg_value
(
    double *x,
    INT     n
) ;

void cg_grad
(
    double  *g,
    double  *x,
    INT      n
) ;

double cg_valgrad
(
    double  *g,
    double  *x,
    INT      n
) ;

/* global variables */
    integer CUTEst_nvar;        /* number of variables */
    integer CUTEst_ncon;        /* number of constraints */

/* main program */
    int MAINENTRY( void ) {

        /* wall clock: */
/*      struct timeval tv ;
        int sec, usec ;
        double walltime ; */
        char *fname = "OUTSDIF.d"; /* CUTEst data file */
        integer funit = 42;        /* FORTRAN unit number for OUTSDIF.d */
        integer io_buffer = 11;    /* FORTRAN unit for internal i/o */
        integer iout = 6;          /* FORTRAN unit number for error output */
        integer ierr;              /* Exit flag from OPEN and CLOSE */
        integer status;            /* Exit flag from CUTEst tools */
        double  grad_tol = 1.e-6; /* required gradient tolerance */

        VarTypes vtypes;

        integer    ncon_dummy ;
        doublereal *x, *bl, *bu ;
        char       *pname, *vnames ;
        logical     efirst = FALSE_, lfirst = FALSE_, nvfrst = FALSE_, grad;
        logical     constrained = FALSE_;

        doublereal  calls[7], cpu[2];
        integer     nlin = 0, nbnds = 0, neq = 0;
        integer     ExitCode;
        int         i, status_cg_descent ;

        FILE *spec ;
        cg_stats Stats ;
        cg_parameter cg_parm ;


        /* Open problem description file OUTSDIF.d */
        ierr = 0;
        FORTRAN_open( &funit, fname, &ierr ) ;
        if( ierr ) {
            printf("Error opening file OUTSDIF.d.\nAborting.\n") ;
            exit(1) ;
        }

        /* Determine problem size */
        CUTEST_cdimen( &status, &funit, &CUTEst_nvar, &CUTEst_ncon) ;
        if (status) {
            printf("** CUTEst error, status = %d, aborting\n", status);
            exit(status);
        }
        /* Determine whether to call constrained or unconstrained tools */
        if( CUTEst_ncon ) constrained = TRUE_;

        /* stop if the problem has constraints */
        if( constrained ) {
           printf (" ** the problem %s has %i constraints\n",
                      pname,  &CUTEst_ncon ) ;
           printf ("    cg_descent is for unconstrained optimization\n") ;
            abort ( ) ;
        }

        /* Seems to be needed for some Solaris C compilers */
        ncon_dummy = CUTEst_ncon + 1;

        /* Reserve memory for variables, bounds, and multipliers */
        /* and call appropriate initialization routine for CUTEst */
        MALLOC( x,  CUTEst_nvar, doublereal ) ;
        MALLOC( bl, CUTEst_nvar, doublereal ) ;
        MALLOC( bu, CUTEst_nvar, doublereal ) ;
        CUTEST_usetup( &status, &funit, &iout, &io_buffer, &CUTEst_nvar,
                       x, bl, bu ) ;
        if (status) {
            printf("** CUTEst error, status = %d, aborting\n", status);
            exit(status);
        }

        /* Get problem name */
        MALLOC( pname,  FSTRING_LEN+1, char );
        CUTEST_probname( &status, pname ) ;
        if (status) {
            printf("** CUTEst error, status = %d, aborting\n", status);
            exit(status);
        }

        /* Make sure to null-terminate problem name */
        pname[FSTRING_LEN] = '\0';
        i = FSTRING_LEN - 1;
        while(i-- > 0 && pname[i] == ' ') {
          pname[i] = '\0';
        }

        printf (" ** the problem is %s\n", pname ) ;

        /* MALLOC(vnames, CUTEst_nvar*FSTRING_LEN, char);
           CUTEST_unames( &status, &CUTEst_nvar, pname, vnames);
           if( status ) {
              printf("** CUTEst error, status = %d, aborting\n", status);
              exit(status);
           }
           FREE(vnames) ;
        */

        /* Set any parameter values here */
        cg_default (&cg_parm) ;

        /* Read input parameters from CG_DESCENT.SPC. See cg_user.h
           for defaults */

        spec = fopen ("CG_DESCENT.SPC", "r") ;

        /*        fscanf( spec, "%lf\n", &grad_tol ) ;*/
        if( fscanf( spec, "%lf\n", &grad_tol )== 1) {
          printf(" grad_tol = %g", grad_tol);
        } else {
          printf("failed to read grad_tol.\n");
         }

        fscanf( spec, "%i\n", &cg_parm.PrintFinal ) ;
        fscanf( spec, "%i\n", &cg_parm.PrintLevel ) ;
        fscanf( spec, "%i\n", &cg_parm.PrintParms ) ;
        fscanf( spec, "%i\n", &cg_parm.LBFGS ) ;
        fscanf( spec, "%i\n", &cg_parm.memory ) ;
        fscanf( spec, "%i\n", &cg_parm.SubCheck ) ;
        fscanf( spec, "%i\n", &cg_parm.SubSkip ) ;
        fscanf( spec, "%lf\n", &cg_parm.eta0 ) ;
        fscanf( spec, "%lf\n", &cg_parm.eta1 ) ;
        fscanf( spec, "%lf\n", &cg_parm.eta2 ) ;
        fscanf( spec, "%i\n", &cg_parm.AWolfe ) ;
        fscanf( spec, "%lf\n", &cg_parm.AWolfeFac ) ;
        fscanf( spec, "%lf\n", &cg_parm.Qdecay ) ;
        fscanf( spec, "%i\n", &cg_parm.nslow ) ;
        fscanf( spec, "%i\n", &cg_parm.StopRule ) ;
        fscanf( spec, "%lf\n", &cg_parm.StopFac ) ;
        fscanf( spec, "%i\n", &cg_parm.PertRule ) ;
        fscanf( spec, "%lf\n", &cg_parm.eps ) ;
        fscanf( spec, "%lf\n", &cg_parm.egrow ) ;
        fscanf( spec, "%i\n", &cg_parm.QuadStep ) ;
        fscanf( spec, "%lf\n", &cg_parm.QuadCutOff ) ;
        fscanf( spec, "%lf\n", &cg_parm.QuadSafe ) ;
        fscanf( spec, "%i\n", &cg_parm.UseCubic ) ;
        fscanf( spec, "%lf\n", &cg_parm.CubicCutOff ) ;
        fscanf( spec, "%lf\n", &cg_parm.SmallCost ) ;
        fscanf( spec, "%i\n", &cg_parm.debug ) ;
        fscanf( spec, "%lf\n", &cg_parm.debugtol ) ;
        fscanf( spec, "%lf\n", &cg_parm.step ) ;
        fscanf( spec, "%li\n", &cg_parm.maxit ) ;
        fscanf( spec, "%i\n", &cg_parm.ntries ) ;
        fscanf( spec, "%lf\n", &cg_parm.ExpandSafe ) ;
        fscanf( spec, "%lf\n", &cg_parm.SecantAmp ) ;
        fscanf( spec, "%lf\n", &cg_parm.RhoGrow ) ;
        fscanf( spec, "%i\n", &cg_parm.neps ) ;
        fscanf( spec, "%i\n", &cg_parm.nshrink ) ;
        fscanf( spec, "%i\n", &cg_parm.nline ) ;
        fscanf( spec, "%lf\n", &cg_parm.restart_fac ) ;
        fscanf( spec, "%lf\n", &cg_parm.feps ) ;
        fscanf( spec, "%lf\n", &cg_parm.nan_rho ) ;
        fscanf( spec, "%lf\n", &cg_parm.nan_decay ) ;
        fscanf( spec, "%lf\n", &cg_parm.delta ) ;
        fscanf( spec, "%lf\n", &cg_parm.sigma ) ;
        fscanf( spec, "%lf\n", &cg_parm.gamma ) ;
        fscanf( spec, "%lf\n", &cg_parm.rho ) ;
        fscanf( spec, "%lf\n", &cg_parm.psi0 ) ;
        fscanf( spec, "%lf\n", &cg_parm.psi_lo ) ;
        fscanf( spec, "%lf\n", &cg_parm.psi_hi ) ;
        fscanf( spec, "%lf\n", &cg_parm.psi1 ) ;
        fscanf( spec, "%lf\n", &cg_parm.psi2 ) ;
        fscanf( spec, "%i\n", &cg_parm.AdaptiveBeta ) ;
        fscanf( spec, "%lf\n", &cg_parm.BetaLower ) ;
        fscanf( spec, "%lf\n", &cg_parm.theta ) ;
        fscanf( spec, "%lf\n", &cg_parm.qeps ) ;
        fscanf( spec, "%lf\n", &cg_parm.qrule ) ;
        fscanf( spec, "%i\n", &cg_parm.qrestart ) ;
        fclose( spec ) ;

/*      cg_parm.debug = TRUE ;*/
/*      cg_parm.PrintLevel = 3 ;*/

        /* Call the optimizer */

        /* wall clock: */
/*      gettimeofday (&tv, NULL) ;
        sec = tv.tv_sec ;
        usec = tv.tv_usec ; */
        status_cg_descent  = cg_descent (x, CUTEst_nvar, &Stats, &cg_parm,
                             grad_tol, cg_value, cg_grad, cg_valgrad, NULL ) ;
/*      gettimeofday (&tv, NULL) ;
        walltime = tv.tv_sec - sec + (double) (tv.tv_usec - usec) /1.e6 ;*/

        ExitCode = 0;

        /* Get CUTEst statistics */
        CUTEST_creport( &status, calls, cpu) ;
        if (status) {
          printf("** CUTEst error, status = %d, aborting\n", status);
          exit(status);
        }

        /* print statistics if so desired */
        printf ("%10s %6i %7li %7li %7li %5i %16.7f %16.7f %9.3f\n",
            pname, CUTEst_nvar, Stats.iter, Stats.nfunc, Stats.ngrad,
            status_cg_descent, Stats.gnorm, Stats.f, cpu [1]) ;
/*          status, Stats.gnorm, Stats.f, walltime) ;*/
/*
        printf("\n\n *********************** CUTEst statistics ************************\n\n") ;
        printf(" Code used               : cg_descent\n") ;
        printf(" Problem                 : %-s\n", pname) ;
        printf(" # variables             = %-10d\n", CUTEst_nvar) ;
        printf(" # bound constraints     = %-10d\n", vtypes.nbnds) ;
        printf(" # objective functions   = %-15.7g\n", calls[0]) ;
        printf(" # objective gradients   = %-15.7g\n", calls[1]) ;
        printf(" # objective Hessians    = %-15.7g\n", calls[2]) ;
        printf(" # Hessian-vector prdct  = %-15.7g\n", calls[3]) ;
        printf(" Exit code               = %-10d\n", ExitCode) ;
        printf(" Final f                 = %-15.7g\n",dummy) ;
        printf(" Set up time             = %-10.2f seconds\n", cpu[0]) ;
        printf(" Solve time              = %-10.2f seconds\n", cpu[1]) ;
        printf(" ******************************************************************\n\n") ;
*/

        ierr = 0;
        FORTRAN_close( &funit, &ierr ) ;
        if( ierr ) {
            printf( "Error closing %s on unit %d.\n", fname, funit ) ;
            printf( "Trying not to abort.\n" ) ;
        }

        /* Free workspace */
        FREE( pname ) ;
        FREE( x ) ; FREE( bl ) ; FREE( bu ) ;

        CUTEST_uterminate( &status ) ;

        return 0;
    }

#ifdef __cplusplus
}    /* Closing brace for  extern "C"  block */
#endif
double cg_value
(
    double *x,
    INT     n
)
{
    double f ;
    integer status;

    CUTEST_ufn( &status, &CUTEst_nvar, x, &f) ;
    if (status) {
        printf("** CUTEst error, status = %d, aborting\n", status);
        exit(status);
    }

    return (f) ;
}

void cg_grad
(
    double  *g,
    double  *x,
    INT      n
)
{
    integer status;
    CUTEST_ugr( &status, &CUTEst_nvar, x, g) ;
    if (status) {
        printf("** CUTEst error, status = %d, aborting\n", status);
        exit(status);
    }
}

double cg_valgrad
(
    double  *g,
    double  *x,
    INT      n
)
{
    logical grad ;
    double f ;
    integer status;
    grad = 1 ;
    CUTEST_uofg( &status, &CUTEst_nvar, x, &f, g, &grad ) ;
    if (status) {
        printf("** CUTEst error, status = %d, aborting\n", status);
        exit(status);
    }
    return (f) ;
}
