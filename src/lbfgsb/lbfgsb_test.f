C     ( Last modified on 4 Jan 2013 at 15:40:00 )

C  Dummy SETULB for testing lbfgsb_main interface to CUTEst
C  Nick Gould,  4th January 2013

      SUBROUTINE SETULB( n, m, X, XL, XU, NBD, f, G, factr, pgtol, WA, 
     *                   IWA, task, iprint, csave, LSAVE, ISAVE, DSAVE )
      INTEGER n, m, iprint, ISAVE( 44 )
      DOUBLE PRECISION f, factr, pgtol, DSAVE( 29 )
      CHARACTER ( LEN = 60 ) :: TASK, CSAVE
      LOGICAL LSAVE( 4 )
      DOUBLE PRECISION NBD( n ), IWA( 3 * n ), X( n ), XL( n ), XU( n ), 
     *    G( n ), WA( 2 * m * n + 4 * n + 12 * m * m + 12 * m )

      IF ( TASK( 1: 5 ) .EQ. 'START' ) THEN
        TASK( 1: 5 ) = 'FG   '
      ELSE
        TASK( 1: 5 ) = 'CONV '
      END IF
      RETURN
      END
