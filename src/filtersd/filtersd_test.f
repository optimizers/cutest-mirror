C     ( Last modified on 30 Jan 2013 at 11:40:00 )

C  Dummy FILTERSD for testing filtersd_main interface to CUTEst
C  Nick Gould,  30th January 2013

      SUBROUTINE filterSD(n, m, x, al, f, fmin, cstype, bl, bu, ws, lws, 
     *  v, nv, maxa, maxla, maxu, maxiu, kmax, maxg, rho, htol, rgtol, 
     *  maxit, iprint, nout, ifail)
      IMPLICIT DOUBLE PRECISION ( a-h, o-z )
      INTEGER :: lws( * )
      DIMENSION :: x( * ), al( * ), bl( * ), bu( * ), ws( * ), v( * )
      CHARACTER :: cstype( * )
      COMMON / defaultc / ainfty, ubd, mlp, mxf
      COMMON / statsc / dnorm, h, hJt, hJ, ipeq, k, itn, nft, ngt
      COMMON / infoc / rgnorm, vstep, iter, npv, ngr, ninf

      last1 = maxu + 1  
      ncx1 = last1 + 2 * maxa
      CALL FUNCTIONS( n, m, x, f, ws( ncx1 ),ws, lws )
      CALL GRADIENTS( n, m, x, ws( last1 ), ws, lws )
      h = 1.0D0
      ubd = 1.0D0
      rgnorm = 1.0D0
      k = n
      itn = 0
      nft = 1
      ngt = 2
      ifail = 5
      RETURN
      END
