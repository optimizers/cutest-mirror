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
        char        fgets_status ;
        char        line[1024];

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

        /* Get problem name (this works under gfortran, but not all compilers*/
        /*
        MALLOC( pname,  FSTRING_LEN+1, char );
        CUTEST_pname( &status, &funit, pname ) ;
        if (status) {
            printf("** CUTEst error, status = %d, aborting\n", status);
            exit(status);
        }
        */
        /* Make sure to null-terminate problem name */
        /*
        pname[FSTRING_LEN] = '\0';
        i = FSTRING_LEN - 1;
        while(i-- > 0 && pname[i] == ' ') {
          pname[i] = '\0';
        }
        */
        /* printf (" ** the problem is %s\n", pname ) ; */

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
                      pname,  CUTEst_ncon ) ;
           printf ("    cg_descent is for unconstrained optimization\n") ;
           exit( -1 ) ;
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

        /*printf ("Problem: %s (n = %i)\n", pname, CUTEst_nvar ) ; */

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

        if ( fgets(line, sizeof line, spec) != NULL ) {
          if( sscanf( line, "%lf\n", &grad_tol ) == 1) {
            /*  printf( " grad_tol = %g\n", grad_tol); */
          } else {
            printf( " failed to read grad_tol\n");
          }
        } else {
          printf( " skipping read grad_tol, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
          if( sscanf( line, "%i\n", &cg_parm.PrintFinal ) == 1) {
            /*  printf( " PrintFinal = %i\n", cg_parm.PrintFinal ) ; */
          } else {
            printf( " failed to read PrintFinal\n");
          }
        } else {
          printf( " skipping read PrintFinal, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%i\n", &cg_parm.PrintLevel ) == 1) {
             /*  printf( " PrintLevel = %i\n", cg_parm.PrintLevel ) ; */
           } else {
             printf( " failed to read PrintLevel\n");
           }
        } else {
         printf( " skipping read PrintLevel, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%i\n", &cg_parm.PrintParms ) == 1) {
             /*  printf( " PrintParms = %i\n", cg_parm.PrintParms ) ; */
           } else {
             printf( " failed to read PrintParms\n");
           }
        } else {
         printf( " skipping read PrintParms, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%i\n", &cg_parm.LBFGS ) == 1) {
             /*  printf( " LBFGS = %i\n", cg_parm.LBFGS ) ; */
           } else {
             printf( " failed to read LBFGS\n");
           }
        } else {
         printf( " skipping read LBFGS, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%i\n", &cg_parm.memory ) == 1) {
             /*  printf( " memory = %i\n", cg_parm.memory ) ; */
           } else {
             printf( " failed to read memory\n");
           }
        } else {
         printf( " skipping read memory, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%i\n", &cg_parm.SubCheck ) == 1) {
             /*  printf( " SubCheck = %i\n", cg_parm.SubCheck ) ; */
           } else {
             printf( " failed to read SubCheck\n");
           }
        } else {
         printf( " skipping read SubCheck, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%i\n", &cg_parm.SubSkip ) == 1) {
             /*  printf( " SubSkip = %i\n", cg_parm.SubSkip ) ; */
           } else {
             printf( " failed to read SubSkip\n");
           }
        } else {
         printf( " skipping read SubSkip, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.eta0 ) == 1) {
             /*  printf( " eta0 = %g\n", cg_parm.eta0 ) ; */
           } else {
             printf( " failed to read eta0\n");
           }
        } else {
         printf( " skipping read eta0, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.eta1 ) == 1) {
             /*  printf( " eta1 = %g\n", cg_parm.eta1 ) ; */
           } else {
             printf( " failed to read eta1\n");
           }
        } else {
         printf( " skipping read eta1, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.eta2 ) == 1) {
             /*  printf( " eta2 = %g\n", cg_parm.eta2 ) ; */
           } else {
             printf( " failed to read eta2\n");
           }
        } else {
         printf( " skipping read eta2, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%i\n", &cg_parm.AWolfe ) == 1) {
             /*  printf( " AWolfe = %i\n", cg_parm.AWolfe ) ; */
           } else {
             printf( " failed to read AWolfe\n");
           }
        } else {
         printf( " skipping read AWolfe, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.AWolfeFac ) == 1) {
             /*  printf( " AWolfeFac = %g\n", cg_parm.AWolfeFac ) ; */
           } else {
             printf( " failed to read AWolfeFac\n");
           }
        } else {
         printf( " skipping read AWolfeFac, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.Qdecay ) == 1) {
             /*  printf( " Qdecay = %g\n", cg_parm.Qdecay ) ; */
           } else {
             printf( " failed to read Qdecay\n");
           }
        } else {
         printf( " skipping read Qdecay, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%i\n", &cg_parm.nslow ) == 1) {
             /*  printf( " nslow = %i\n", cg_parm.nslow ) ; */
           } else {
             printf( " failed to read nslow\n");
           }
        } else {
         printf( " skipping read nslow, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%i\n", &cg_parm.StopRule ) == 1) {
             /*  printf( " StopRule = %i\n", cg_parm.StopRule ) ; */
           } else {
             printf( " failed to read StopRule\n");
           }
        } else {
         printf( " skipping read StopRule, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.StopFac ) == 1) {
             /*  printf( " StopFac = %g\n", cg_parm.StopFac ) ; */
           } else {
             printf( " failed to read StopFac\n");
           }
        } else {
         printf( " skipping read StopFac, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%i\n", &cg_parm.PertRule ) == 1) {
             /*  printf( " PertRule = %i\n", cg_parm.PertRule ) ; */
           } else {
             printf( " failed to read PertRule\n");
           }
        } else {
         printf( " skipping read PertRule, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.eps ) == 1) {
             /*  printf( " eps = %g\n", cg_parm.eps ) ; */
           } else {
             printf( " failed to read eps\n");
           }
        } else {
         printf( " skipping read eps, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.egrow ) == 1) {
             /*  printf( " egrow = %g\n", cg_parm.egrow ) ; */
           } else {
             printf( " failed to read egrow\n");
           }
        } else {
         printf( " skipping read egrow, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%i\n", &cg_parm.QuadStep ) == 1) {
             /*  printf( " QuadStep = %i\n", cg_parm.QuadStep ) ; */
           } else {
             printf( " failed to read QuadStep\n");
           }
        } else {
         printf( " skipping read QuadStep, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.QuadCutOff ) == 1) {
             /*  printf( " QuadCutOff = %g\n", cg_parm.QuadCutOff ) ; */
           } else {
             printf( " failed to read QuadCutOff\n");
           }
        } else {
         printf( " skipping read QuadCutOff, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.QuadSafe ) == 1) {
             /*  printf( " QuadSafe = %g\n", cg_parm.QuadSafe ) ; */
           } else {
             printf( " failed to read QuadSafe\n");
           }
        } else {
         printf( " skipping read QuadSafe, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%i\n", &cg_parm.UseCubic ) == 1) {
             /*  printf( " UseCubic = %i\n", cg_parm.UseCubic ) ; */
           } else {
             printf( " failed to read UseCubic\n");
           }
        } else {
         printf( " skipping read UseCubic, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.CubicCutOff ) == 1) {
             /*  printf( " CubicCutOff = %g\n", cg_parm.CubicCutOff ) ; */
           } else {
             printf( " failed to read CubicCutOff\n");
           }
        } else {
         printf( " skipping read CubicCutOff, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.SmallCost ) == 1) {
             /*  printf( " SmallCost = %g\n", cg_parm.SmallCost ) ; */
           } else {
             printf( " failed to read SmallCost\n");
           }
        } else {
         printf( " skipping read SmallCost, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%i\n", &cg_parm.debug ) == 1) {
             /*  printf( " debug = %i\n", cg_parm.debug ) ; */
           } else {
             printf( " failed to read debug\n");
           }
        } else {
         printf( " skipping read debug, status = %c\n", fgets_status );
        }
       if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.debugtol ) == 1) {
             /*  printf( " debugtol = %g\n", cg_parm.debugtol ) ; */
           } else {
             printf( " failed to read debugtol\n");
           }
        } else {
         printf( " skipping read debugtol, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.step ) == 1) {
             /*  printf( " step = %g\n", cg_parm.step ) ; */
           } else {
             printf( " failed to read step\n");
           }
        } else {
         printf( " skipping read step, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%li\n", &cg_parm.maxit ) == 1) {
             /*  printf( " maxit = %i\n", cg_parm.maxit ) ; */
           } else {
             printf( " failed to read maxit\n");
           }
        } else {
         printf( " skipping read maxit, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%i\n", &cg_parm.ntries ) == 1) {
             /*  printf( " ntries = %i\n", cg_parm.ntries ) ; */
           } else {
             printf( " failed to read ntries\n");
           }
        } else {
         printf( " skipping read ntries, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.ExpandSafe ) == 1) {
             /*  printf( " ExpandSafe = %g\n", cg_parm.ExpandSafe ) ; */
           } else {
             printf( " failed to read ExpandSafe\n");
           }
        } else {
         printf( " skipping read ExpandSafe, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.SecantAmp ) == 1) {
             /*  printf( " SecantAmp = %g\n", cg_parm.SecantAmp ) ; */
           } else {
             printf( " failed to read SecantAmp\n");
           }
        } else {
         printf( " skipping read SecantAmp, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.RhoGrow ) == 1) {
             /*  printf( " RhoGrow = %g\n", cg_parm.RhoGrow ) ; */
           } else {
             printf( " failed to read RhoGrow\n");
           }
        } else {
         printf( " skipping read RhoGrow, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%i\n", &cg_parm.neps ) == 1) {
             /*  printf( " neps = %i\n", cg_parm.neps ) ; */
           } else {
             printf( " failed to read neps\n");
           }
        } else {
         printf( " skipping read neps, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%i\n", &cg_parm.nshrink ) == 1) {
             /*  printf( " nshrink = %i\n", cg_parm.nshrink ) ; */
           } else {
             printf( " failed to read nshrink\n");
           }
        } else {
         printf( " skipping read nshrink, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%i\n", &cg_parm.nline ) == 1) {
             /*  printf( " nline = %i\n", cg_parm.nline ) ; */
           } else {
             printf( " failed to read nline\n");
           }
        } else {
         printf( " skipping read nline, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.restart_fac ) == 1) {
             /*  printf( " restart_fac = %g\n", cg_parm.restart_fac ) ; */
           } else {
             printf( " failed to read restart_fac\n");
           }
        } else {
         printf( " skipping read restart_fac, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.feps ) == 1) {
             /*  printf( " feps = %g\n", cg_parm.feps ) ; */
           } else {
             printf( " failed to read feps\n");
           }
        } else {
         printf( " skipping read feps, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.nan_rho ) == 1) {
             /*  printf( " nan_rho = %g\n", cg_parm.nan_rho ) ; */
           } else {
             printf( " failed to read nan_rho\n");
           }
        } else {
         printf( " skipping read nan_rho, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.nan_decay ) == 1) {
             /*  printf( " nan_decay = %g\n", cg_parm.nan_decay ) ; */
           } else {
             printf( " failed to read nan_decay\n");
           }
        } else {
         printf( " skipping read nan_decay, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.delta ) == 1) {
             /*  printf( " delta = %g\n", cg_parm.delta ) ; */
           } else {
             printf( " failed to read delta\n");
           }
        } else {
         printf( " skipping read delta, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.sigma ) == 1) {
             /*  printf( " sigma = %g\n", cg_parm.sigma ) ; */
           } else {
             printf( " failed to read sigma\n");
           }
        } else {
         printf( " skipping read sigma, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.gamma ) == 1) {
             /*  printf( " gamma = %g\n", cg_parm.gamma ) ; */
           } else {
             printf( " failed to read gamma\n");
           }
        } else {
         printf( " skipping read gamma, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.rho ) == 1) {
             /*  printf( " rho = %g\n", cg_parm.rho ) ; */
           } else {
             printf( " failed to read rho\n");
           }
        } else {
         printf( " skipping read rho, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.psi0 ) == 1) {
             /*  printf( " psi0 = %g\n", cg_parm.psi0 ) ; */
           } else {
             printf( " failed to read psi0\n");
           }
        } else {
         printf( " skipping read psi0, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.psi_lo ) == 1) {
             /*  printf( " psi_lo = %g\n", cg_parm.psi_lo ) ; */
           } else {
             printf( " failed to read psi_lo\n");
           }
        } else {
         printf( " skipping read psi_lo, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.psi_hi ) == 1) {
             /*  printf( " psi_hi = %g\n", cg_parm.psi_hi ) ; */
           } else {
             printf( " failed to read psi_hi\n");
           }
        } else {
         printf( " skipping read psi_hi, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.psi1 ) == 1) {
             /*  printf( " psi1 = %g\n", cg_parm.psi1 ) ; */
           } else {
             printf( " failed to read psi1\n");
           }
        } else {
         printf( " skipping read psi1, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.psi2 ) == 1) {
             /*  printf( " psi2 = %g\n", cg_parm.psi2 ) ; */
           } else {
             printf( " failed to read psi2\n");
           }
        } else {
         printf( " skipping read psi2, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%i\n", &cg_parm.AdaptiveBeta ) == 1) {
             /*  printf( " AdaptiveBeta = %i\n", cg_parm.AdaptiveBeta ) ; */
           } else {
             printf( " failed to read AdaptiveBeta\n");
           }
        } else {
         printf( " skipping read AdaptiveBeta, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.BetaLower ) == 1) {
             /*  printf( " BetaLower = %g\n", cg_parm.BetaLower ) ; */
           } else {
             printf( " failed to read BetaLower\n");
           }
        } else {
         printf( " skipping read BetaLower, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.theta ) == 1) {
             /*  printf( " theta = %g\n", cg_parm.theta ) ; */
           } else {
             printf( " failed to read theta\n");
           }
        } else {
         printf( " skipping read theta, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.qeps ) == 1) {
             /*  printf( " qeps = %g\n", cg_parm.qeps ) ; */
           } else {
             printf( " failed to read qeps\n");
           }
        } else {
         printf( " skipping read qeps, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%lf\n", &cg_parm.qrule ) == 1) {
             /*  printf( " qrule = %g\n", cg_parm.qrule ) ; */
           } else {
             printf( " failed to read qrule\n");
           }
        } else {
         printf( " skipping read qrule, status = %c\n", fgets_status );
        }
        if ( fgets(line, sizeof line, spec) != NULL ) {
           if( sscanf( line, "%i\n", &cg_parm.qrestart ) == 1) {
             /*  printf( " qrestart = %i\n", cg_parm.qrestart ) ; */
           } else {
             printf( " failed to read qrestart\n");
           }
        } else {
         printf( " skipping read qrestart, status = %c\n", fgets_status );
        }

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
        /*
        printf ("%10s %6i %7li %7li %7li %5i %16.7f %16.7f %9.3f\n",
            pname, CUTEst_nvar, Stats.iter, Stats.nfunc, Stats.ngrad,
            status_cg_descent, Stats.gnorm, Stats.f, cpu [1]) ;
        */
/*          status, Stats.gnorm, Stats.f, walltime) ;*/
        printf(" *********************** CUTEst statistics ************************\n") ;
        printf(" Code used               : cg_descent\n") ;
        printf(" Problem                 : %-s\n", pname) ;
        printf(" # variables             = %-10d\n", CUTEst_nvar) ;
        printf(" # bound constraints     = %-10d\n", vtypes.nbnds) ;
        printf(" # iterations            = %li\n", Stats.iter) ;
        printf(" # objective functions   = %-15.7g\n", calls[0]) ;
        printf(" # objective gradients   = %-15.7g\n", calls[1]) ;
        printf(" # objective Hessians    = %-15.7g\n", calls[2]) ;
        printf(" # Hessian-vector prdct  = %-15.7g\n", calls[3]) ;
        printf(" Exit code               = %-10d\n", ExitCode) ;
        printf(" Final f                 = %-15.7g\n",Stats.f) ;
        printf(" Final ||g||             = %-15.7g\n",Stats.gnorm) ;
        printf(" Set up time             = %-10.2f seconds\n", cpu[0]) ;
        printf(" Solve time              = %-10.2f seconds\n", cpu[1]) ;
        printf(" ******************************************************************\n") ;

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
