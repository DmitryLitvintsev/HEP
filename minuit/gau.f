      FUNCTION GAU(X)
C*
C* LAMC 2PI SPECTRUM 
C*
      COMMON/PAWPAR/PAR(11)
      REAL M0,M,M1,M2,NSIGM,LEFT,RIGHT

      AMPL1 =     PAR(1)
      SIGMA1=     PAR(2)
      M1    =     PAR(3) 

      A0   =     PAR(4)
      A1   =     PAR(5)
      A2   =     PAR(6)
      A3   =     PAR(7)


C      AMPL1  = -63.244
C      SIGMA1 = 8.0323
C      M1     = 1.8612 
C      A0     = 0.20637E+06
C      A1     =  -0.28017E+06
C      A2     = 0.10579E+06
C      A3     = 62749.
C
      AMPL1  = -18.021
      SIGMA1 = 7.8919
      M1     =  1.8618
      A0     =  69061.
      A1     = -0.10434E+06
      A2     = 65482.
      A3     = -8032.8
C
      BIN=5.
      ATHR=1.321+0.1396

      
      DM=X-ATHR
      

      GAU=SQRT(0.5*(ABS(DM)+DM))*
     &       (A0+A1*DM+A2*DM**2+A3*DM**3)*5.E-3+
     &       AMPL1*BIN*BWG(X,0.0,SIGMA1,M1)
      END
C*
C* CONVOLUTION
C*
      FUNCTION BWG(X,GAMMA,SIGMA,M0)
      REAL M0,M,NSIGM,LEFT,RIGHT,GAMMA,X,SIGMA,VAR,SUM
      EPSILON=GAMMA/2/SIGMA
      M=1000.*X
      M0=1000.*M0
c     M=X
      NSIGM = (M-M0)/SIGMA  

      IF (GAMMA.LT.0.0001) THEN
         BWG=0.3989423/SIGMA*EXP(-0.5*NSIGM**2) 
        GOTO 200
      ENDIF

      NPOINTS=1000
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





