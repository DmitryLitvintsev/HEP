      FUNCTION FIT1(X)
C*
C* LAMC 2PI SPECTRUM 
C*
      COMMON/PAWPAR/PAR(8)
      REAL M0,M,M1,M2,NSIGM,LEFT,RIGHT

      AMPL1 =     PAR(1)
      SIGMA1=     PAR(2)
      GAMMA1=     PAR(3)
      M1    =     PAR(4) 

      A0    =     PAR(5)
      A1    =     PAR(6)
      A2    =     PAR(7)
      A3    =     PAR(8)

      BIN   = 2.E-3
      ATHR  = 1.321+0.1396
      DM    = X-ATHR
      

      FIT1=SQRT(0.5*(ABS(DM)+DM))*
     &       (A0+A1*DM+A2*DM**2+A3*DM**3)+
     &       AMPL1*BIN*BWG(X,GAMMA1,SIGMA1,M1)
      END
C*
C* CONVOLUTION
C*
C
C This is handmade convolution
C
      FUNCTION BWG(X,GAMMA,SIGMA,M0)
      REAL M0,M,NSIGM,LEFT,RIGHT,GAMMA,X,SIGMA,VAR,SUM
      REAL EPSILON,STEP
      INTEGER NPOINTS

      EPSILON = 0.5*GAMMA/SIGMA
      M       = X
      NSIGM = (M-M0)/SIGMA  

      IF (GAMMA.LT.0.0001) THEN
        FOLD1=0.3989423/SIGMA*DEXP(-0.5*NSIGM**2) 
        GOTO 201
      ENDIF

      NPOINTS = 100
      LEFT    = (-5.+NSIGM)/EPSILON
      RIGHT   = (5.+NSIGM)/EPSILON
      STEP    = (RIGHT-LEFT)/FLOAT(NPOINTS)
      SUM=0.
      DO 101 I=1,NPOINTS
       VAR=LEFT+(FLOAT(I)-0.5)*STEP
 101   SUM=SUM+DEXP(-0.5*(NSIGM-EPSILON*VAR)**2)/(VAR**2+1.)
       BWG=SUM*STEP*2./SIGMA/15.7496
 201   CONTINUE
      RETURN
      END







