/* ============================================
 * CUTEst interface for GNU Scientific Library (GSL)
 *
 * J. Hogg 2015
 * Based on GENC interface
 * ============================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl_blas.h>
#include <gsl_matrix.h>
#include <gsl_multifit_nlin.h>
#include <gsl_vector.h>

#define GENCMA

#ifdef __cplusplus
extern "C" {   /* To prevent C++ compilers from mangling symbols */
#endif

#include "cutest.h"

#define GENC    genc
#define GENSPC  genspc
#define GETINFO getinfo

double GENC( doublereal );
void GENSPC( integer, char * );
void GETINFO( integer, integer, doublereal *, doublereal *,
              doublereal *, doublereal *, logical *, logical *,
              VarTypes * );

integer CUTEst_nvar;        /* number of variables */
integer CUTEst_ncon;        /* number of constraints */
integer CUTEst_nnzj;        /* number of nonzeros in Jacobian */
integer CUTEst_nnzh;        /* number of nonzeros in upper triangular
                                  part of the Hessian of the Lagrangian */

struct param_data {
   int m;
   int n;
};

int eval_fn( const gsl_vector *x, void *params, gsl_vector *f ) {
   struct param_data *data = (struct param_data*) params;
   int status;
   double obj, i;
   const double *xptr = gsl_vector_const_ptr(x, 0);
   double *fptr = gsl_vector_ptr(f, 0);
   CUTEST_cfn( &status, &data->n, &data->m, xptr, &obj, fptr);
   return GSL_SUCCESS;
}
int eval_jacobian( const gsl_vector *x, void *params, gsl_matrix *J ) {
   struct param_data *data = (struct param_data*) params;
   double *y = (double *) malloc(data->m*sizeof(double));
   double *g = (double *) malloc(data->n*sizeof(double));
   bool grlagf;
   bool jtrans;
   const double *xptr;
   double *Jptr;
   int status, ldJ;
   grlagf = false;
   jtrans = true; /* GSL uses row major */
   xptr = gsl_vector_const_ptr(x, 0);
   Jptr = gsl_matrix_ptr(J, 0, 0);
   ldJ = J->tda;
   CUTEST_cgr( &status, &data->n, &data->m, xptr, y, &grlagf, g, &jtrans,
         &ldJ, &data->m, Jptr );
   free(y); free(g); /* Values ignored */
   return GSL_SUCCESS;
}
int eval_fn_jacobian( const gsl_vector *x, void *params, gsl_vector *f,
      gsl_matrix *J ) {
   struct param_data *data = (struct param_data*) params;
   double *y = (double *) malloc(data->m*sizeof(double));
   double *g = (double *) malloc(data->n*sizeof(double));
   bool grlagf;
   bool jtrans;
   const double *xptr;
   double *fptr, *Jptr;
   int status, ldJ;
   double obj;
   grlagf = false;
   jtrans = true; /* GSL uses row major */
   xptr = gsl_vector_const_ptr(x, 0);
   Jptr = gsl_matrix_ptr(J, 0, 0);
   fptr = gsl_vector_ptr(f, 0);
   CUTEST_cfn( &status, &data->n, &data->m, xptr, &obj, fptr);
   ldJ = J->tda;
   CUTEST_cgr( &status, &data->n, &data->m, xptr, y, &grlagf, g, &jtrans,
         &ldJ, &data->m, Jptr );
   free(y); free(g); /* Values ignored */
   return GSL_SUCCESS;
}

int MAINENTRY( void ){

    char *fname = "OUTSDIF.d"; /* CUTEst data file */
    integer funit = 42;        /* FORTRAN unit number for OUTSDIF.d */
    integer iout = 6;          /* FORTRAN unit number for error output */
    integer io_buffer = 11;    /* FORTRAN unit internal input/output */
    integer ierr;              /* Exit flag from OPEN and CLOSE */
    integer status;            /* Exit flag from CUTEst tools */
    integer gsl_status;        /* Exit flag from GSL */
    integer grad_status;       /* Exit flag from computing gradient */
    
    FILE *indr, *summary;

    VarTypes vtypes;

    doublereal *x, *bl, *bu, *dummy1, *dummy2;
    doublereal *v = NULL, *cl = NULL, *cu = NULL;
    logical *equatn = NULL, *linear = NULL;
    char *pname, *vnames, *gnames, *cptr;
    char **Vnames; /* vnames and gnames as arrays of strings */
    logical grad;
    integer e_order = 1, l_order = 0, v_order = 0;
    logical constrained = FALSE_;

    doublereal calls[7], cpu[2];
    integer nlin = 0, nbnds = 0, neq = 0;
    integer ExitCode;
    int i, j;

    int maxiter, iter;
    double epsabs, epsrel;
    struct param_data params;
    gsl_vector_view xview;
    gsl_vector *gradient;
    gsl_multifit_function_fdf fdf;
    gsl_multifit_fdfsolver *gsl;
    int info;
    double fnval;

    int write_summary, summary_size;
    int stopping_test;
    double normJf, normf;
    double tol_gradient, tol_func;
    char summary_file[20];

    /* Open problem description file OUTSDIF.d */
    ierr = 0;
    FORTRAN_open( &funit, fname, &ierr );
    if ( ierr )
    {
        printf("Error opening file OUTSDIF.d.\nAborting.\n");
        exit(1);
    }

    /* Determine problem size */
    CUTEST_cdimen( &status, &funit, &CUTEst_nvar, &CUTEst_ncon );

    if ( status )
    {
        printf("** CUTEst error, status = %d, aborting\n", status);
        exit(status);
    }

    /* Determine whether to call constrained or unconstrained tools */
    if ( CUTEst_ncon ) constrained = TRUE_;

    /* Reserve memory for variables, bounds, and multipliers */
    /* and call appropriate initialization routine for CUTEst */
    MALLOC( x,      CUTEst_nvar, doublereal );
    MALLOC( bl,     CUTEst_nvar, doublereal );
    MALLOC( bu,     CUTEst_nvar, doublereal );
    if ( constrained )
    {
        MALLOC( equatn, CUTEst_ncon, logical    );
        MALLOC( linear, CUTEst_ncon, logical    );
        MALLOC( v,      CUTEst_ncon, doublereal );
        MALLOC( cl,     CUTEst_ncon, doublereal );
        MALLOC( cu,     CUTEst_ncon, doublereal );
        CUTEST_csetup( &status, &funit, &iout, &io_buffer,
                       &CUTEst_nvar, &CUTEst_ncon, x, bl, bu,
                       v, cl, cu, equatn, linear,
                       &e_order, &l_order, &v_order );
        /*        printf("CUTEst_nvar = %d\n", CUTEst_nvar);
        printf("CUTEst_ncon = %d\n", CUTEst_ncon);
        printf("x = ");
        for (i = 0; i < CUTEst_nvar ; i++)
            printf("%g ", x[i]);
        printf("\n");
        printf("bl = ");
        for (i = 0; i < CUTEst_nvar ; i++)
            printf("%g ", bl[i]);
        printf("\n");
        printf("bu = ");
        for (i = 0; i < CUTEst_nvar ; i++)
            printf("%g ", bu[i]);
        printf("\n");
        printf("v = ");
        for (i = 0; i < CUTEst_ncon ; i++)
            printf("%g ", v[i]);
        printf("\n");
        printf("cl = ");
        for (i = 0; i < CUTEst_ncon ; i++)
            printf("%g ", cl[i]);
        printf("\n");
        printf("cu = ");
        for (i = 0; i < CUTEst_ncon ; i++)
            printf("%g ", cu[i]);
        printf("\n");
        printf("equatn = ");
        for (i = 0; i < CUTEst_ncon ; i++)
            printf("%d ", equatn[i]);
        printf("\n");
        printf("linear = ");
        for (i = 0; i < CUTEst_ncon ; i++)
            printf("%d ", linear[i]);
            printf("\n"); */
    }
    else
        CUTEST_usetup( &status, &funit, &iout, &io_buffer,
                       &CUTEst_nvar, x, bl, bu );

    if ( status )
    {
        printf("** CUTEst error, status = %d, aborting\n", status);
        exit(status);
    }

    /* Get problem, variables and constraints names */
    MALLOC(pname, FSTRING_LEN + 1, char);
    MALLOC(vnames, CUTEst_nvar * FSTRING_LEN, char);        /* For Fortran */
    MALLOC(Vnames, CUTEst_nvar, char *);               /* Array of strings */
    for (i = 0; i < CUTEst_nvar; i++)
        MALLOC(Vnames[i], FSTRING_LEN + 1, char);

    if ( constrained )
    {
        MALLOC(gnames, CUTEst_ncon * FSTRING_LEN, char);   /* For Fortran */
        CUTEST_cnames( &status, &CUTEst_nvar, &CUTEst_ncon,
                       pname, vnames, gnames );
    }
    else
    {
        CUTEST_unames( &status, &CUTEst_nvar, pname, vnames );
    }

    if ( status )
    {
        printf("** CUTEst error, status = %d, aborting\n", status);
        exit(status);
    }

    /* Make sure to null-terminate problem name */
    pname[FSTRING_LEN] = '\0';

    /* Transfer variables and constraint names into arrays of
     * null-terminated strings.
     * If you know of a simpler way to do this portably, let me know!
     */
    for (i = 0; i < CUTEst_nvar; i++)
    {
        cptr = vnames + i * FSTRING_LEN;
        for (j = 0; j < FSTRING_LEN; j++)
        {
            Vnames[i][j] = *cptr;
            cptr++;
        }
        Vnames[i][FSTRING_LEN] = '\0';
    }

    /* Fortran strings no longer needed */
    FREE(vnames);
    if (constrained) FREE(gnames);

    /* Obtain basic info on problem */
    GETINFO(CUTEst_nvar, CUTEst_ncon, bl, bu, cl, cu,
            equatn, linear, &vtypes);

    printf("Initial point:\n");
    for (i = 0; i < CUTEst_nvar; i++)
        printf("  %s %e\n", Vnames[i], x[i]);

    /* Initialize optimizer */
    gsl = gsl_multifit_fdfsolver_alloc(
             gsl_multifit_fdfsolver_lmsder, CUTEst_ncon, CUTEst_nvar
             );
    /* User defined parameter data */
    params.m = CUTEst_ncon;
    params.n = CUTEst_nvar;
    /* function information */
    fdf.f = eval_fn;
    fdf.df = eval_jacobian;
    fdf.fdf = eval_fn_jacobian;
    fdf.n = CUTEst_ncon;
    fdf.p = CUTEst_nvar;
    fdf.params = &params;
    /* Set up optimizer with function and initial guess */
    xview = gsl_vector_view_array(x, CUTEst_nvar);
    gsl_multifit_fdfsolver_set(gsl, &fdf, &xview.vector);

    /* Convergence parameters */
    /* Stop if |dx_i| < epsabs + epsrel |x_i| for all i. */
    /*    maxiter = 1000;
    epsabs = 1e-12;
    epsrel = 1e-12;*/
    
    /* Read parameters in from GSL.SPC */
    indr = fopen("GSL.SPC","r"); 
    if ( indr == NULL ) {
      fprintf(stderr, "Error: can't find GLS.SPC\n");
      exit(1);
    }
    rewind( indr );
    
    ierr = fscanf( indr, "%i%*[^\n]\n", &maxiter);
    if (ierr != 1) {
      printf("Error: failed to read max iterations from GSL.SPC; using default (1000). \n");
      maxiter = 1000;
    }
    ierr = fscanf( indr, "%lf%*[^\n]\n", &epsabs);
    if (ierr != 1) {
      printf("ierr = %i",ierr);
      printf("Error: failed to read abs tolerance from GSL.SPC; using default (1e-6). \n");
      epsabs = 1e-6;
    }
    ierr = fscanf( indr, "%lf%*[^\n]\n", &epsrel);
    if (ierr != 1) {
      printf("Error: failed to read abs tolerance from GSL.SPC; using default (1e-6). \n");
      epsrel = 1e-6;
    }
    ierr = fscanf( indr, "%i%*[^\n]\n", &stopping_test);
    if (ierr != 1) {
      printf("Error: failed to read stopping_test from GSL.SPC; using default (1). \n");
      write_summary = 1;
    }
    ierr = fscanf( indr, "%i%*[^\n]\n", &write_summary);
    if (ierr != 1) {
      printf("Error: failed to read write_summary from GSL.SPC; using default (0). \n");
      write_summary = 0;
    }
    ierr = fscanf( indr, "%s", summary_file);
    if (ierr != 1) {
      printf("Error: failed to read summary_file from GSL.SPC; using default (GSL.SPC). \n");
      strncpy(summary_file,"GSL.SPC",20);
    }
    
    fclose(indr);

    if ( write_summary == 1) {
      printf("Writing summary to %s \n", summary_file);
      summary = fopen(summary_file,"a"); /* append summary, create if doesn't exist */
      if ( summary == NULL ) {
	fprintf(stderr, "Error: can't open %s for writing\n",summary_file);
	exit(1);
      }
      summary_size = ftell(summary);
      if (summary_size == 0) { 
	/* write where opened if file was empty */
	  fprintf( summary, "%10s", pname );
	} 
      else {
	/* otherwise start with a new line */
	fprintf( summary, "\n%10s", pname );
      }
    }
    
    /* Call the optimizer */
    /*ExitCode = gsl_multifit_fdfsolver_driver(gsl, maxiter, epsabs, epsrel);
    switch(ExitCode) {
        case GSL_SUCCESS:
            printf("Converged sucessfully (info=%d)\n", info);
            break;
        case GSL_EMAXITER:
            printf("Reached maximum number of iterations\n");
            break;
        default:
            printf("Unknown error %d (info=%d)\n", ExitCode, info);
            break;
    }*/
    iter = 0;

    /* first, get the initial values... */
    gradient = gsl_vector_alloc(fdf.p);
    normf = gsl_blas_dnrm2(gsl->f);
    fnval = 0.5*normf*normf;
    printf("0.5 ||f0||^2 = %e\n", fnval);
    grad_status = gsl_multifit_gradient(gsl->J, gsl->f, gradient);
    normJf = gsl_blas_dnrm2(gradient);
    printf("||J0'f0||/||f0|| = %e\n", normJf/normf);
    
    if (stopping_test == 1) {

      /* set tol_gradient = max(epsabs, ||J^Tf||/||f|| * epsrel) */
      tol_gradient = (normJf/normf) * epsrel;
      if (tol_gradient < epsabs){
	tol_gradient = epsabs;
      }
      
      /* set tol_func = max(epsabs, ||f|| * epsrel ) */
      tol_func = normf * epsrel;
      if (tol_func < epsabs) { 
	tol_func = epsabs;
      }
            
    }

    /* now, to the main loop... */
    do {
       printf("===== Iteration %d =====\n", iter);
       iter++;
       gsl_status = gsl_multifit_fdfsolver_iterate(gsl);
       printf("itr status =%s\n", gsl_strerror(gsl_status));
       if(gsl_status) break;

       /* check for convergence */       
       normf = gsl_blas_dnrm2(gsl->f);
       fnval = 0.5*normf*normf;
       printf("0.5 ||f||^2 = %e\n", fnval);

       
       grad_status = gsl_multifit_gradient(gsl->J, gsl->f, gradient);
       normJf = gsl_blas_dnrm2(gradient);
       printf("||J'f||/||f|| = %e\n", normJf/normf);
       
	  
       if (grad_status) break;
       if (stopping_test == 1) {
	 /* test on the scaled gradient */
	 if ( (normJf/normf) > tol_gradient) {
	   gsl_status = GSL_CONTINUE;
	 }
	 else {
	   gsl_status = GSL_SUCCESS;
	   printf("Scaled gradient test successful\n");
	   break;
	 }
	 
	 /* test on the function value */
	 if (normf > tol_func ){
	   gsl_status = GSL_CONTINUE;
	 }
	 else {
	   gsl_status = GSL_SUCCESS;
	   printf("Norm of function test successful\n");
	   break;
	 }

       }
       else if (stopping_test == 2) {
	 /* Stop if |dx_i| < epsabs + epsrel |x_i| for all i. */
	 gsl_status = gsl_multifit_test_delta(gsl->dx, gsl->x, epsabs, epsrel);
       } 
       else if (stopping_test == 3) { 
	 /* stop if || g ||_inf < epsabs */
	 gsl_status = gsl_multifit_test_gradient(gradient,epsabs);
       } 
       else {
	 printf("Error: unsupported stopping test \n");
	 gsl_status = GSL_EINVAL; /* made up error code... */
       }
       printf("test status = %s\n", gsl_strerror(gsl_status));
    }
    while(gsl_status==GSL_CONTINUE && iter < maxiter);

    /* Calculate final value */
    /*    normf = gsl_blas_dnrm2(gsl->f);
	  fnval = 0.5*normf*normf;*/

    /* Get CUTEst statistics */
    CUTEST_creport( &status, calls, cpu );

    if ( status )
    {
        printf("** CUTEst error, status = %d, aborting\n", status);
        exit(status);
    }

    /* write summary if required */
    if ( write_summary == 1) {
      fprintf( 
	      summary, "%6i %6i %6i %6i %23.15e %23.15e",
	      (int)CUTEst_nvar,
	      (int)CUTEst_ncon,
	      gsl_status, 
	      iter, 
	      fnval,
	      normJf
	);
    }

    /* print stats to screen */

    printf("===== Results =====\n");
    printf("Variables at optimum:\n");
    for (i = 0; i < CUTEst_nvar; i++)
        printf("  %s %e\n", Vnames[i], gsl_vector_get(gsl->x, i));
    printf("Final function value = %e\n", fnval);

    printf("\n\n ************************ CUTEst statistics ************************\n\n");
    printf(" Code used               : GSL\n");
    printf(" Problem                 : %-s\n", pname);
    printf(" # variables             = %-10d\n", (int)CUTEst_nvar);
    printf(" # constraints           = %-10d\n", (int)CUTEst_ncon);
    printf(" # linear constraints    = %-10d\n", vtypes.nlin);
    printf(" # equality constraints  = %-10d\n", vtypes.neq);
    printf(" # inequality constraints= %-10d\n", vtypes.nineq);
    printf(" # bound constraints     = %-10d\n", vtypes.nbnds);
    printf(" # objective functions   = %-15.7g\n", calls[0]);
    printf(" # objective gradients   = %-15.7g\n", calls[1]);
    printf(" # objective Hessians    = %-15.7g\n", calls[2]);
    printf(" # Hessian-vector prdct  = %-15.7g\n", calls[3]);
    if (constrained) printf(" # constraints functions = %-15.7g\n", calls[4]);
    if (constrained) printf(" # constraints gradients = %-15.7g\n", calls[5]);
    if (constrained) printf(" # constraints Hessians  = %-15.7g\n", calls[6]);
    printf(" Exit code               = %-10d\n", (int)ExitCode);
    printf(" Final f                 = %-15.7e\n", fnval);
    printf(" Set up time             = %-10.2f seconds\n", cpu[0]);
    printf(" Solve time              = %-10.2f seconds\n", cpu[1]);
    printf(" ******************************************************************\n\n");

    ierr = 0;
    FORTRAN_close( &funit, &ierr );
    if ( ierr )
    {
        printf( "Error closing %s on unit %d.\n", fname, (int)funit );
        printf( "Trying not to abort.\n" );
    }

    /* Free workspace */
    FREE( pname );
    FREE( x ); FREE( bl ); FREE( bu );
    FREE( v ); FREE( cl ); FREE( cu );
    FREE( equatn );
    FREE( linear );
    for (i = 0; i < CUTEst_nvar; i++) FREE(Vnames[i]);
    FREE(Vnames);

    if ( constrained )
      CUTEST_cterminate( &status );
    else
      CUTEST_uterminate( &status );

    return 0;

}

void getinfo( integer n, integer m, doublereal *bl, doublereal *bu,
              doublereal *cl, doublereal *cu, logical *equatn,
              logical *linear, VarTypes *vartypes )
{

    int i;

    vartypes->nlin = 0; vartypes->neq = 0; vartypes->nbnds = 0;
    vartypes->nrange = 0;
    vartypes->nlower = 0; vartypes->nupper = 0; vartypes->nineq = 0;
    vartypes->nineq_lin = 0; vartypes->nineq_nlin = 0;
    vartypes->neq_lin = 0; vartypes->neq_nlin = 0;

    for ( i = 0; i < n; i++ )
        if ( bl[i] > -CUTE_INF || bu[i] < CUTE_INF ) vartypes->nbnds++;
    for ( i = 0; i < m; i++ )
    {
        if ( linear[i] ) vartypes->nlin++;
        if ( equatn[i] )
        {
            vartypes->neq++;
            if ( linear[i] )
                vartypes->neq_lin++;
            else
                vartypes->neq_nlin++;
        }
        else
        {
            vartypes->nineq++;
            if ( cl[i] > -CUTE_INF )
            {
                if ( cu[i] < CUTE_INF )
                    vartypes->nrange++;
                else
                {
                    vartypes->nlower++;
                }
            }
            else
            {
                if ( cu[i] < CUTE_INF )
                {
                    vartypes->nupper++;
                }
            }
            if ( linear[i] )
            {
                vartypes->nineq_lin++;
            }
            else
            {
                vartypes->nineq_nlin++;
            }
        }
    }
    return;
}


#ifdef __cplusplus
}    /* Closing brace for  extern "C"  block */
#endif