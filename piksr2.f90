!**************************************************************************
      SUBROUTINE PIKSR2 (N, ARR, BRR)
!**************************************************************************
!
!     Try the simple insertion sort.
!
!**************************************************************************
      REAL ARR(N), A
      INTEGER BRR(N), B
      
      DO J = 2, N
        A = ARR(J)
        B = BRR(J)
        DO I = J-1,1,-1
          IF(ARR(I).LE.A) GOTO 10
          ARR(I+1)=ARR(I)
          BRR(I+1)=BRR(I)
        END DO
        I = 0
  10    ARR(I+1)=A
        BRR(I+1)=B
      END DO
      
      RETURN
      END
