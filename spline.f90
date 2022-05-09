!***********************************************************************
      SUBROUTINE SPLINE (N,X,Y,FDP)      
!***********************************************************************
!-----THIS SUBROUTINE COMPUTES THE SECOND DERIVATIVES NEEDED 
!-----IN CUBIC SPLINE INTERPOLATION.  THE INPUT DATA ARE:    
!-----N = NUMBER OF DATA POINTS          
!-----X = ARRAY CONTAINING THE VALUES OF THE INDEPENDENT VARIABLE      
!-----    (ASSUMED TO BE IN ASCENDING ORDER)       
!-----Y = ARRAY CONTAINING THE VALUES OF THE FUNCTION AT THE 
!-----    DATA POINTS GIVEN IN THE X ARRAY         
!-----THE OUTPUT IS THE ARRAY FDP WHICH CONTAINS THE SECOND  
!-----DERIVATIVES OF THE INTERPOLATING CUBIC SPLINE.         
      DIMENSION X(N),Y(N),A(N),B(N),C(N),R(N),FDP(N)  
!-----COMPUTE THE COEFFICIENTS AND THE RHS OF THE EQUATIONS. 
!-----THIS ROUTINE USES THE CANTILEVER CONDITION.  THE PARAMETER       
!-----ALAMDA (LAMBDA) IS SET TO 1. BUT THIS CAN BE USER-MODIFIED.      
!-----A,B,C ARE THE THREE DIAGONALS OF THE TRIDIAGONAL SYSTEM;         
!-----R IS THE RIGHT HAND SIDE.  THESE ARE NOW ASSEMBLED.    
      ALAMDA = 1.    
      NM2 = N - 2    
      NM1 = N - 1    
      C(1) = X(2) - X(1)       
      DO 1 I=2,NM1   
      C(I) = X(I+1) - X(I)     
      A(I) = C(I-1)  
      B(I) = 2.*(A(I) + C(I))  
      R(I) = 6.*((Y(I+1) - Y(I))/C(I) - (Y(I) - Y(I-1))/C(I-1))        
    1 CONTINUE       
      B(2) = B(2) + ALAMDA * C(1)        
      B(NM1) = B(NM1) + ALAMDA * C(NM1)  
!-----AT THIS POINT WE COULD CALL A TRIDIAGONAL SOLVER SUBROUTINE      
!-----BUT THE NOTATION IS CLUMSY SO WE WILL SOLVE DIRECTLY.  THE       
!-----NEXT SECTION SOLVES THE SYSTEM WE HAVE JUST SET UP.    
      DO 2 I=3,NM1   
      T = A(I)/B(I-1)          
      B(I) = B(I) - T * C(I-1) 
      R(I) = R(I) - T * R(I-1) 
    2 CONTINUE       
      FDP(NM1) = R(NM1)/B(NM1) 
      DO 3 I=2,NM2   
      NMI = N - I    
      FDP(NMI) = (R(NMI) - C(NMI)*FDP(NMI+1))/B(NMI)         
    3 CONTINUE       
      FDP(1) = ALAMDA * FDP(2) 
      FDP(N) = ALAMDA * FDP(NM1)         
!-----WE NOW HAVE THE DESIRED DERIVATIVES SO WE RETURN TO THE          
!-----MAIN PROGRAM.  
      RETURN         
      END  

!***********************************************************************
      SUBROUTINE SPEVAL (N,X,Y,FDP,XX,F) 
!***********************************************************************
!-----THIS SUBROUTINE EVALUATES THE CUBIC SPLINE GIVEN       
!-----THE DERIVATIVE COMPUTED BY SUBROUTINE SPLINE.          
!-----THE INPUT PARAMETERS N,X,Y,FDP HAVE THE SAME 
!-----MEANING AS IN SPLINE.    
!-----XX = VALUE OF INDEPENDENT VARIABLE FOR WHICH 
!-----     AN INTERPOLATED VALUE IS REQUESTED      
!-----F =  THE INTERPOLATED RESULT       
      DIMENSION X(N),Y(N),FDP(N)      
!-----THE FIRST JOB IS TO FIND THE PROPER INTERVAL.          
      NM1 = N - 1    
      DO 1 I=1,NM1   
      IF (XX.LE.X(I+1)) GO TO 10         
    1 CONTINUE       
!-----NOW EVALUATE THE CUBIC   
   10 DXM = XX - X(I)          
      DXP = X(I+1) - XX        
      DEL = X(I+1) - X(I)      
      F = FDP(I)*DXP*(DXP*DXP/DEL - DEL)/6.     &
          + FDP(I+1)*DXM*(DXM*DXM/DEL - DEL)/6. &     
          + Y(I)*DXP/DEL + Y(I+1)*DXM/DEL
      RETURN        
      END 
