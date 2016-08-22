C
C unbinned likelihood fit to fix Xi(1530) and Xi(1862) only
C
      PROGRAM MAIN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /PAWC/ LPAW(300000)
      INTEGER MXDATA,NDATA
      PARAMETER(MXDATA=1000000)
      DOUBLE PRECISION XDATA,EXDATA
      COMMON /USER/ XDATA(MXDATA),EXDATA(MXDATA),NDATA
      INTEGER ISTAT,NUMBER,RUN,EVENT
      REAL RNDM,DUMMY
      EXTERNAL FCN

      CALL HLIMIT(300000)
      
      CALL HROPEN(20,'blah','tst.hist','N',1024,istat)
      CALL HBOOK1(1,'xicst',700.0,1.4,2.8,0.0)
      CALL HBOOK1(2,'LH',2000,-1000.0,1000.0,0.0)


      NDATA=1
      OPEN(UNIT=80,FILE='data_points_deltaR_rs_refit.dat')
 50   READ(80,*,END=500)NUMBER,RUN,EVENT,
     &     XDATA(NDATA),EXDATA(NDATA)
      NDATA=NDATA+1
      GOTO 50
 500  NDATA=NDATA-1
      GOTO 999
 800  FORMAT(1X,'Total data pts read:  ',i10)
 999  WRITE(*,800)NDATA
      CLOSE(80)

C*
C*    Add test gaussian
C*
      XMASS = 1.862
      SIGMA = 0.008 



C      I = 1 
C      DO WHILE(I.LE.1000)
C        XM      = (1.D0-2.D0*RNDM(DUMMY))*SIGMA*10.D0
C        PRB     = RNDM(DUMMY)
C        Y       = XM / SIGMA
C        GAU     = DEXP(-0.5*Y**2)
C        IF (PRB .LE. GAU) THEN
C          XDATA(NDATA+I) = XMASS + XM
C          I = I + 1 
C        ENDIF
C      ENDDO
C      NDATA = NDATA + 1000
C      I = 1 
      DO I=1,NDATA
        CALL HF1(1,FLOAT(XDATA(I)),1.0)
      ENDDO

C*
C*   Open data files
C*
      open(unit=5,file='xicst2.data',status='old')
      open(unit=21,file='fcn.dat',status='unknown',form='formatted')
      open(unit=22,file='func.dat',status='unknown',form='formatted')
C      open(unit=6,file='fit.out',status='unknown',form='formatted')

      CALL MINTIO(5,6,7)
      CALL MINUIT(FCN,0)

C      CALL HRPUT(0,'tst.hist','N')
      CALL HROUT (0,ICYCLE,' ')
      CALL HREND('blah')
      
      CLOSE(5)
      CLOSE(21)
      CLOSE(22)

      STOP
      END


      SUBROUTINE FCN(NPAR,GIN,F,U,IFLAG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(*),GIN(*)
      INTEGER MXDATA,NDATA
      PARAMETER(MXDATA=1000000)
      DOUBLE PRECISION XDATA,EXDATA
      COMMON /USER/ XDATA(MXDATA),EXDATA(MXDATA),NDATA
      DOUBLE PRECISION FUN,SUM,SQRT,DEXP,BWG,C1,C2,HIGHM
      DOUBLE PRECISION GCM,GCS,GCG

      HIGHM = 2.2D0
      C1    = HIGHM - (1.321D0+0.1396D0)

      IF (IFLAG.NE.1) GOTO 100
      
 100  SUM=U(1)+U(9)+DSQRT(C1)*(U(5)*2.D0/3.D0*C1+
     &  U(6)*2.D0/5.D0*C1**2+
     &  U(7)*2.D0/7.D0*C1**3+
     &  U(8)*2.D0/9.D0*C1**4)
C
C Gaussian Constraints
C
      GCM   = (U(11)-1.862D0)/0.002D0
      GCS   = (U(10)-0.008D0)/0.0012D0
      DO I=1,NDATA
        IF (XDATA(I).LE.HIGHM) THEN
          FUN= BWG(XDATA(I),U)
          IF (FUN.LE.0) THEN
            SUM = SUM + 1000.0
          ELSE
            SUM=SUM-DLOG(FUN)
          ENDIF
        ENDIF
      END DO
      SUM = SUM + 0.5*GCS*GCS + 0.5*GCM*GCM
      F=SUM*2.D0
      IF (IFLAG.EQ.3) THEN
        WRITE(21,1000  )F
 1000   FORMAT( G20.10 )
        
        
C        SUM=U(1)+U(9)+DSQRT(C1)*(U(5)*2.D0/3.D0*C1+
C     &       U(6)*2.D0/5.D0*C1**2+
C     &       U(7)*2.D0/7.D0*C1**3+
C     &       U(8)*2.D0/9.D0*C1**4)
C        
C        DO IDATA=1,NDATA
C          FUN= BWG(XDATA(IDATA),U)
C          WRITE(22,*) FUN/SUM
C        END DO
C        
C
C        RINT = 0.D0
C        DO I=-1000,1000
C          U(9)  = DFLOAT(I)
C          SUM=U(1)+U(9)+DSQRT(C1)*(U(5)*2.D0/3.D0*C1+
C     &         U(6)*2.D0/5.D0*C1**2+
C     &         U(7)*2.D0/7.D0*C1**3+
C     &         U(8)*2.D0/9.D0*C1**4)
C          DO IDATA=1,NDATA
C            IF (XDATA(IDATA).LE.HIGHM) THEN
C              FUN= BWG(XDATA(IDATA),U)
C              IF (FUN.LE.0) THEN
C                SUM = SUM + 1000.0
C              ELSE
C                SUM=SUM-DLOG(FUN)
C              ENDIF
C            ENDIF
C          END DO
C          SUM = SUM + 0.5*GCS*GCS + 0.5*GCM*GCM
C          RF=SUM*2.D0
C          RF = RF-F
C          CALL HF1(2,U(9),RF)
C          WRITE(21,1001  )RF
C 1001     FORMAT( G20.10 )
C        ENDDO
      ENDIF
      RETURN
      END
      

 
      FUNCTION GAUSS(OBM,TRM,SIG)
      IMPLICIT NONE 
      DOUBLE PRECISION X,OBM,TRM,SIG,GAUSS,DSQRT,DEXP
      X = (OBM - TRM) / SIG
      GAUSS = 1.D0/DSQRT(2.D0*3.14159265358979312D0)/SIG*
     &       DEXP(-0.5D0*X**2)
      RETURN
      END
C
C homemade function
C
      
      FUNCTION FOLD1(X,GAMMA,SIGMA,M0)
      IMPLICIT NONE 
      DOUBLE PRECISION M0,M,NSIGM,LEFT,RIGHT,GAMMA,X,SIGMA,VAR,SUM
      DOUBLE PRECISION EPSILON,DEXP,FOLD1,STEP
      INTEGER NPOINTS,I
      EPSILON=0.5D0*GAMMA/SIGMA
      M = X
      NSIGM = (M-M0)/SIGMA  

      IF (GAMMA.LT.0.0001D0) THEN
        FOLD1=0.3989423D0/SIGMA*DEXP(-0.5D0*NSIGM**2) 
        GOTO 201
      ENDIF

      NPOINTS=1000
      LEFT  = (-5.D0+NSIGM)/EPSILON
      RIGHT = ( 5.D0+NSIGM)/EPSILON
      STEP  = (RIGHT-LEFT)/DFLOAT(NPOINTS)
      SUM=0.D0
      DO 101 I=1,NPOINTS
        VAR=LEFT+(DFLOAT(I)-0.5D0)*STEP
 101    SUM=SUM+DEXP(-0.5D0*(NSIGM-EPSILON*VAR)**2)/(VAR**2+1.D0)
        FOLD1=SUM*STEP*2.D0/SIGMA/15.7496D0
 201    CONTINUE
      RETURN
      END

C
C CERN based function
C

      FUNCTION FOLD(OBM,GAM,SIG,TRM)
      IMPLICIT NONE 
      COMPLEX*16 Z,WZ,WWERF
      DOUBLE PRECISION X,Y,FOLD,OBM,GAM,SIG,TRM
      X = (OBM - TRM) / (SIG * 1.414213D0)
      Y = GAM / (2.828427D0 * SIG)
      Z = CMPLX(X,Y)
      WZ = WWERF(Z)
      FOLD = DREAL(WZ)
      FOLD = FOLD / (2.50662827D0 * SIG )
      RETURN
      END
 
      FUNCTION BWG(X,U)
      IMPLICIT NONE
      DOUBLE PRECISION X,U
      DOUBLE PRECISION N1,SIG1,GAMMA,M1,A0,A1,A2,A3,DM
      DOUBLE PRECISION N2,SIG2,M2
      DOUBLE PRECISION BWG,FOLD,DSQRT,GAUSS,FOLD1
      DIMENSION U(*)
      


      N1       = U(1)
      SIG1     = U(2)
      GAMMA    = U(3)
      M1       = U(4)

      A0       = U(5)
      A1       = U(6)
      A2       = U(7)
      A3       = U(8)

      N2       = U(9)
      SIG2     = U(10)
      M2       = U(11)


      DM = X - (1.321D0+0.1396D0)

      IF (DM.LT.0) THEN
        BWG = 0.D0
        RETURN
      ENDIF

      BWG=N1*FOLD(X,GAMMA,SIG1,M1)+
     &    N2*GAUSS(X,M2,SIG2)+
     &     DSQRT(DM)*
     &     (A0+A1*DM+A2*DM**2+A3*DM**3)
      
      RETURN
      END
      
      FUNCTION PS(X,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IF (X.GT.0.D0 .AND. N.GT.1) THEN
         S1 = DFLOAT(N) * DLOG(X)
         S2 = 0.D0
         DO 10 I = 2,N
10       S2 = S2 + DLOG(DFLOAT(I))
         PS = S1 - S2 - X
      ENDIF
      IF (X.EQ.0.D0             ) PS= 1.D0
      IF (X.GT.0.D0 .AND. N.EQ.0) PS=-1.D0 * X
      IF (X.GT.0.D0 .AND. N.EQ.1) PS=DLOG(X) - X
      IF (X.LT.0                ) THEN
          PRINT *,'******* FUNCTION LESS THAN 0 ********'
          PRINT *,'******* PS IS SET TO -100    ********'
          PS = -100.D0
      ENDIF
      RETURN
      END

