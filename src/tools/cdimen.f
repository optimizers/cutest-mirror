! ( Last modified on 10 Sepc 2004 at 16:35:38 )
!  Correction: 10/Sep/2004: undeclared integer variables declared
      SUBROUTINE CDIMEN( INPUT, N, M )
      INTEGER :: INPUT, N, M

!  Compute the basic array dimensions for the constrained optimization tools.

!  Nick Gould, for CGT productions,
!  26th August, 1999.

      INTEGER :: IALGOR, I, J, NG, NELNUM, IDUMMY
      INTEGER :: NG1, NEL1, NSLACK, NOBJGR, IEND
      INTEGER :: IARRAY( 10 )
      CHARACTER ( LEN = 8 ) :: PNAME

!  Input the problem dimensions.

      REWIND INPUT
      READ( INPUT, 1001 ) N, NG, NELNUM

!  Input the problem type.

      READ( INPUT, 1000 ) IALGOR, PNAME
      IF ( IALGOR < 2 ) THEN
         M = 0
         GO TO 100
      END IF

!  Set useful integer values.

      NG1 = NG + 1
      NEL1 = NELNUM + 1

!  Print out problem data. input the number of variables, groups,
!  elements and the identity of the objective function group.

      IF ( IALGOR == 2 ) READ( INPUT, 1002 ) NSLACK, NOBJGR

!  Input the starting addresses of the elements in each group,
!  of the parameters used for each group and
!  of the nonzeros of the linear element in each group.

      READ( INPUT, 1010 ) ( IDUMMY, I = 1, NG1 )
      READ( INPUT, 1010 ) ( IDUMMY, I = 1, NG1 )
      READ( INPUT, 1010 ) ( IDUMMY, I = 1, NG1 )

!  Input the starting addresses of the variables and parameters
!  in each element.

      READ( INPUT, 1010 ) ( IDUMMY, I = 1, NEL1 )
      READ( INPUT, 1010 ) ( IDUMMY, I = 1, NEL1 )

!  Input the group type of each group

      READ( INPUT, 1010 ) ( IDUMMY, I = 1, NG )

!  Count the number of constraint groups

      M = 0
      DO 20 I = 1, NG, 10
         IEND = MIN( I + 9, NG )
         READ( INPUT, 1010 ) ( IARRAY( J - I + 1 ), J = I, IEND )
         DO 10 J = I, IEND
           IF ( IARRAY( J - I + 1 ) /= 1 ) M = M + 1
   10    CONTINUE
   20 CONTINUE

  100 CONTINUE
      REWIND INPUT
      RETURN

!  Non-executable statements.

 1000 FORMAT( I2, A8 )
 1001 FORMAT( 3I8 )
 1002 FORMAT( 2I8 )
 1010 FORMAT( ( 10I8 ) )

!  End of CDIMEN.

      END
