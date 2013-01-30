A grand unified Matlab gateway for the CUTEst tools.
This interface brings together the unconstrained, constrained,
dense and sparse versions of the CUTEst tools.

In order to unify the tools and be able to use the same Matlab commands on
both constrained and unconstrained problems, the tool names in this
interface differ from the those in the old Fortran gateway routine.

--------------------------------------------------------------------------
Matlab Tool  CUTEst library function(s)  Purpose
--------------------------------------------------------------------------
dims         cdimen                      Obtain problem dimensions
setup        usetup / csetup             Setup problem data structure
obj          uofg / cofg                 Evaluate objective function value
                                          and its gradient if requested

objcons      cfn                         Evaluate objective and constraints

cons         ccfg / ccifg                Evaluate constraint bodies
                                          and their gradients if requested.
                                         Evaluate a single constraint value
                                          and its gradient if requested

scons        ccfsg / ccifsg              Evaluate constraint bodies and
                                          Jacobian in sparse format.
                                         Evaluate a single constraint value
                                          and its gradient as a sparse vector

lagjac       cgr                         Evaluate Jacobian and gradient of
                                          either objective or Lagrangian

slagjac      csgr                        Evaluate Jacobian in sparse format
                                          and gradient of either objective or
                                          Lagrangian as a sparse vector

Jprod        cjprod                      Evaluate the matrix-vector product
                                          between the Jacobian and a vector

Jtprod       cjprod                      Evaluate the matrix-vector product
                                          between the transpose Jacobian and
                                          a vector

hess         udh / cdh                   Evaluate the Hessian matrix of the
                                          Lagrangian, or of the objective if
                                          the problem is unconstrained

ihess        udh / cidh                  Evaluate the Hessian matrix of the
                                          i-th problem function (i=0 is the
                                          objective function), or of the
                                          objective if problem is unconstrained

hprod        uprod / cprod               Evaluate the matrix-vector product
                                          between the Hessian of the
                                          Lagrangian
                                          (or the objective if unconstrained)
                                          and a vector

gradhess      ugrdh / cgrdh              Evaluate the gradient of either the
                                          objective or the Lagrangian, the
                                          Jacobian (or its transpose) and the
                                          Hessian of the Lagrangian in dense
                                          format

sphess       ush / csh                   Evaluate the Hessian matrix of the
                                          Lagrangian, or of the objective if
                                          the problem is unconstrained, in
                                          sparse format

isphess      ush / cish                  Evaluate the Hessian matrix of the
                                          i-th problem function (i=0 is the
                                          objective function), or of the
                                          objective if problem is
                                          unconstrained, in sparse format
--------------------------------------------------------------------------

* To create the Matlab/CUTEst interface, use the cutest2matlab command or,
  for more control, the runcutest command with the  -p matlab option. 
  This will create a binary file mcutest.mexglx (32bit Linux), 
  mcutest.mexa64 (64bit Linux) mcutest.mexmaci (32bit OSX) or 
  mcutest.mexmaci64 (64bit OSX).

* Both this binary and the Matlab .m files in $CUTEST/src/matlab must be 
  placed on the Matlab search path. 

* See the man page for runcutest for more details of other options, and the 
  Matlab help files for the tools mentioned above.

CUTEr version:
 D. Orban, Montreal, January 2007
CUTEst version additions:
 Nick Gould, January 2013
