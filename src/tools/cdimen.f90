! ( Last modified on 10 Sepc 2004 at 16:35:38 )
!  Correction: 10/Sep/2004: undeclared integer variables declared
      SUBROUTINE CDIMEN( input, n, m )
      INTEGER :: input, n, m

!  Compute the basic array dimensions for the constrained optimization tools.

!  Nick Gould, for CGT productions,
!  26th August, 1999.

      INTEGER :: ialgor, i, j, ng, nel, idummy
      INTEGER :: ng1, nel1, nslack, nobjgr, iend
      INTEGER :: IARRAY( 10 )
      CHARACTER ( LEN = 8 ) :: PNAME

!  Input the problem dimensions.

      REWIND input
      READ( input, 1001 ) n, ng, nel

!  Input the problem type.

      READ( input, 1000 ) ialgor, PNAME
      IF ( ialgor < 2 ) THEN
         m = 0
         GO TO 100
      END IF

!  Set useful integer values.

      ng1 = ng + 1
      nel1 = nel + 1

!  Print out problem data. input the number of variables, groups,
!  elements and the identity of the objective function group.

      IF ( ialgor == 2 ) READ( input, 1002 ) nslack, nobjgr

!  Input the starting addresses of the elements in each group,
!  of the parameters used for each group and
!  of the nonzeros of the linear element in each group.

      READ( input, 1010 ) ( idummy, i = 1, ng1 )
      READ( input, 1010 ) ( idummy, i = 1, ng1 )
      READ( input, 1010 ) ( idummy, i = 1, ng1 )

!  Input the starting addresses of the variables and parameters
!  in each element.

      READ( input, 1010 ) ( idummy, i = 1, nel1 )
      READ( input, 1010 ) ( idummy, i = 1, nel1 )

!  Input the group type of each group

      READ( input, 1010 ) ( idummy, i = 1, ng )

!  Count the number of constraint groups

      m = 0
      DO 20 i = 1, ng, 10
         iend = MIN( i + 9, ng )
         READ( input, 1010 ) ( IARRAY( j - i + 1 ), j = i, iend )
         DO 10 j = i, iend
           IF ( IARRAY( j - i + 1 ) /= 1 ) m = m + 1
   10    CONTINUE
   20 CONTINUE

  100 CONTINUE
      REWIND input
      RETURN

!  Non-executable statements.

 1000 FORMAT( I2, A8 )
 1001 FORMAT( 3I8 )
 1002 FORMAT( 2I8 )
 1010 FORMAT( ( 10I8 ) )

!  End of CDIMEN.

      END
