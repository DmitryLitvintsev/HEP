      REAL FUNCTION FIT(X)
C*
C* LAMC 2PI SPECTRUM 
C*
      COMMON/PAWPAR/PAR(8)
      REAL M0,M,M1,M2,NSIGM,LEFT,RIGHT

      AMPL1 =     PAR(1)
      SIGMA1=     PAR(2)
      GAMMA1=     PAR(3)
      M1    =     PAR(4) 

      A0   =     PAR(5)
      A1   =     PAR(6)
      A2   =     PAR(7)
      A3   =     PAR(8)

C      A1 = 1000
C      SIGMA1 = 3.2
C      GAMMA1 = 9.1 
C      M1 = 1.531

C      A0 = 0
C      A1 = 0 
C      A2 = 0 
C      A3 = 0 
      BIN=2.
      ATHR=1.321+0.1396
      DM=X-M1
      

      FIT=SQRT(0.5*(ABS(X-ATHR)+X-ATHR))*
     &       (A0+A1*DM+A2*DM**2+A3*DM**3)+
     &       AMPL1*BIN*BWG(X,ABS(GAMMA1),SIGMA1,M1)
      END
C*
C* CONVOLUTION
C*
      FUNCTION BWG(X,GAMMA,SIGMA,M0)
      REAL M0,M,NSIGM,LEFT,RIGHT,GAMMA,X,SIGMA,VAR,SUM
      EPSILON=GAMMA/2/SIGMA
      M=1000.*X
c     M=X
      NSIGM = (M-M0)/SIGMA  

      IF (GAMMA.LT.0.0001) THEN
         BWG=0.3989423/SIGMA*EXP(-0.5*NSIGM**2) 
        GOTO 200
      ENDIF

      NPOINTS=100
      LEFT=(-5+NSIGM)/EPSILON
      RIGHT=(5+NSIGM)/EPSILON
      STEP=(RIGHT-LEFT)/NPOINTS
      SUM=0.
      DO 100 I=1,NPOINTS
       VAR=LEFT+(I-0.5)*STEP
 100   SUM=SUM+EXP(-0.5*(NSIGM-EPSILON*VAR)**2)/(VAR**2+1)
       BWG=SUM*STEP*2./SIGMA/15.7496
 200  CONTINUE
      RETURN
      END





