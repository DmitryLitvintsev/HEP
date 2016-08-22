      FUNCTION FIT2(X)
C*
C* LAMC 2PI SPECTRUM 
C*
      COMMON/PAWPAR/PAR(11)
      REAL M0,M,M1,M2,NSIGM,LEFT,RIGHT

      AMPL1 =     PAR(1)
      SIGMA1=     PAR(2)
      GAMMA1=     PAR(3)
      M1    =     PAR(4) 

      A0   =     PAR(5)
      A1   =     PAR(6)
      A2   =     PAR(7)
      A3   =     PAR(8)

      AMPL2   = PAR(9)
      SIGMA2  = PAR(10)
      M2      = PAR(11)
      
      SIG3    = PAR(12)
      R1      = PAR(13)
      SIG4    = PAR(14)
      R2      = PAR(15)
      
C FCN=  -1162269.     FROM MIGRAD    STATUS=CONVERGED    409 CALLS     2083 TOTAL
C                     EDM=  0.43E-04    STRATEGY= 2      ERROR MATRIX ACCURATE 
C
C  EXT PARAMETER                                   STEP         FIRST   
C  NO.   NAME        VALUE          ERROR          SIZE      DERIVATIVE 
C   1    N1          1887.8        70.193        35.055       0.69117E-06
C   2    SIG1       0.37000E-02     fixed    
C   3    GAMMA      0.91000E-02   constant   
C   4    M1          1.5317       0.29813E-03   0.15688E-03   0.18499    
C   5    A0         0.25956E+06    2137.4        278.75       0.55754E-07
C   6    A1        -0.47539E+06    8470.0        524.66      -0.20860E-06
C   7    A2         0.38163E+06    15379.        845.22      -0.41303E-07
C   8    A3         -67509.        16448.        1251.3       0.12194E-06
C   9    N2          47.154        42.324        21.009      -0.56179E-05
C  10    SIG2       0.39003E-02   0.79988E-04   0.42101E-04    1.1122    
C  11    M2          1.8627       0.19920E-02   0.10146E-02   -4.3592    
C  12    SIG3       0.15000E-02     fixed    
C  13    R1          3.2075        3.0335        1.5324       0.60714E-03
C  14    SIG4       0.96004E-02   0.20000E-03   0.10528E-03   0.79982E-01
C  15    R2          1.9394       0.16004       0.84233E-01  -0.70692E-03

      AMPL1  = 1887.8
      SIGMA1 = 3.7000
      GAMMA1 = 9.1 
      M1     = 1.5317
      A0     = 0.25956E+06
      A1     = -0.47539E+06
      A2     = 0.38163E+06
      A3     = -67509.

      M2      =  1.8627
      SIGMA2  =  3.900
      AMPL2   =  47.154 

      SIG3    = 1.5
C      R1      = 0.
      R1      = 3.2075 
      SIG4    = 9.6
      R2      = 1.9394

C      AMPL1  = 1886.2
C      SIGMA1 = 3.2727
C      GAMMA1 = 9.1 
C      M1     = 1.5317
C      A0     =  0.25978E+06
C      A1     = -0.47745E+06 
C      A2     = 0.38616E+06
C      A3     = -70371.
C
C      M2      = 1.8621 
C      SIGMA2  = 8.0412
C      AMPL2   = 59.179 
C

C      AMPL1  = 676.79
C      SIGMA1 = 3.0437
C      GAMMA1 = 9.1
C      M1     = 1.5321
C      A0     =  86439.
C      A1     = -0.20573E+06
C      A2     = 0.26716E+06 
C      A3     = -0.13549E+06
C
C      M2      = 1.8620 
C      SIGMA2  = 8.0010
C      AMPL2   = 0.22414
C
      BIN=5.
      ATHR=1.321+0.1396

      
      DM=X-ATHR
      

C      FIT2=SQRT(0.5*(ABS(DM)+DM))*
C     &       (A0+A1*DM+A2*DM**2+A3*DM**3)*5.E-3+
C     &       AMPL1*BIN*(R1/(1.+R1)*BWG(X,ABS(GAMMA1),SIGMA1,M1)+
C     &       1.0/(1.+R1)*BWG(X,ABS(GAMMA1),SIG3,M1))+
C     &       AMPL2*BIN*(R2/(1.+R2)*BWG(X,0.0,SIGMA2,M2)+
C     &       1.0/(1.+R2)*BWG(X,0.0,SIG4,M2))

      FIT2=SQRT(0.5*(ABS(DM)+DM))*
     &       (A0+A1*DM+A2*DM**2+A3*DM**3)*5.E-3+
     &       AMPL1*BIN*BWG(X,ABS(GAMMA1),SIG3,M1)+
     &       AMPL2*BIN*(R2/(1.+R2)*BWG(X,0.0,SIGMA2,M2)+
     &       1.0/(1.+R2)*BWG(X,0.0,SIG4,M2))
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





