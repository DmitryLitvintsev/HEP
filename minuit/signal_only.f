C
C this is to run many experiments and check how they look
C
      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER MXDATA,NDATA,NPARM
      INTEGER LPAW
      COMMON /PAWC/ LPAW(300000)
      PARAMETER(MXDATA=1000000)
      DOUBLE PRECISION XDATA,EXDATA
      COMMON /USER/ XDATA(MXDATA),EXDATA(MXDATA),NDATA
      DOUBLE PRECISION MASS,MASS0,SIGMA,GAMMA,PRB,FOLD,FOLD1
      INTEGER ISTAT,NUMBER,RUN,EVENT, NEXP, I, J, IDATA, ICYCLE
      REAL RNDM2,DUMMY,RNDM
      EXTERNAL FCN
      PARAMETER(NPARM=4)
      CHARACTER*10 PNAM(NPARM) 
      INTEGER NPRM
      DOUBLE PRECISION VSTRT, STP, ARGLIS
      DIMENSION NPRM(NPARM),VSTRT(NPARM),STP(NPARM),ARGLIS(10) 
      INTEGER IERFLG
      DATA NPRM /   1   ,    2   ,     3    ,   4    / 
      DATA PNAM / 'N1', 'SIG1', 'GAMMA','M1' /
      DATA VSTRT/  2000.D0  , 3.2D-3  ,   9.1D-3 , 1.5318D0 /
      DATA STP  /   0.001 ,    0.000 ,     0.001  ,   0.001 /



      CALL HLIMIT(300000)
      CALL HROPEN(20,'blah','signal1.hist','N',1024,ISTAT)
      CALL HBOOK1(1,'mass',560,1.4,2.8,0.0)
      
 
      MASS0  = 1.5318D0
      NDATA  = 2000 
      NEXP   = 1000
      SIGMA  = 3.2D-3
      GAMMA  = 9.1D-3
      IDATA  = 0
      I      = 0

      
C*     
C*   Open data files
C*     
C      open(unit=5,file='signal_only.data',status='old')
      open(unit=22,file='parms.dat',
     &     status='unknown',form='formatted')
      open(unit=23,file='residuals.dat',
     &     status='unknown',form='formatted')
      open(unit=6,file='parms.out',
     &     status='unknown',form='formatted')
      
      CALL MNINIT(5,6,7) 

      DO I=1,NPARM
        CALL MNPARM(NPRM(I),PNAM(I),VSTRT(I),STP(I),0.D0,0.D0,IERFLG) 
        IF (IERFLG .NE. 0)  THEN 
          WRITE (6,'(A,I)')  ' UNABLE TO DEFINE PARAMETER NO.',I 
          STOP 
        ENDIF  
      ENDDO

      I = 0 
      
      CALL MNSETI('xicst')
      
 1    CONTINUE
      IDATA = 0
      DO WHILE ( IDATA .LT. NDATA ) 
C        MASS = MASS0 + (1.D0-2.D0*RNDM())*10.D0*5.6D-3
        MASS = RNDM()*3.
        PRB  = RNDM()*FOLD1(MASS0,GAMMA,SIGMA,MASS0) 
        IF (PRB.LE.FOLD1(MASS,GAMMA,SIGMA,MASS0)) THEN
          IDATA = IDATA+1
          XDATA(IDATA) = MASS 
          CALL HF1(1,REAL(MASS),1.0)
        ENDIF
      ENDDO       

      ARGLIS(1) = 1. 
      CALL MNEXCM(FCN, 'CALL FCN', ARGLIS ,1,IERFLG,0) 

      ARGLIS(1) = 2. 
C      CALL MNEXCM(FCN,'FIX', ARGLIS ,1,IERFLG,0) 
      CALL MNEXCM(FCN,'SET PRINT', ARGLIS ,1,IERFLG,0) 
      CALL MNEXCM(FCN,'SET STR', ARGLIS ,1,IERFLG,0) 
      ARGLIS(1) = 0.
      CALL MNEXCM(FCN,'MIGRAD', ARGLIS ,0,IERFLG,0) 
C      CALL MNEXCM(FCN,'MINOS', ARGLIS ,0,IERFLG,0) 
      ARGLIS(1) = 3.
      CALL MNEXCM(FCN,'CALL FCN', ARGLIS , 1,IERFLG,0) 
      CALL MNEXCM(FCN,'STOP ', 0.D0,0,IERFLG,0)

C      CALL MINUIT(FCN,0)
      I=I+1
      IF (I.LT.NEXP) GOTO 1

      
      CLOSE(5)
C     CLOSE(6)
      CLOSE(21)
      CLOSE(22)
      CALL HROUT (0,ICYCLE,' ')
      CALL HREND('blah')
      
 
      STOP
      END
      
      
      SUBROUTINE FCN(NPAR,GIN,F,U,IFLAG)
      IMPLICIT NONE
      INTEGER IFLAG,NPAR
      DOUBLE PRECISION U,GIN,F
      DIMENSION U(*),GIN(*)
      INTEGER MXDATA,NDATA,NPARM,I,NPOINTS
      PARAMETER(MXDATA=1000000)
      DOUBLE PRECISION XDATA,EXDATA
      COMMON /USER/ XDATA(MXDATA),EXDATA(MXDATA),NDATA
      DOUBLE PRECISION FUN,SUM,SQRT,DEXP,BWG,BWG2,C1,C2,HIGHM
      DOUBLE PRECISION GCM,GCS,GCG
      PARAMETER(NPARM=4)
      DOUBLE PRECISION EU
      DIMENSION EU(NPARM)
      DOUBLE PRECISION EMAT,EPLUS,EMINUS,EPARAB,GLOBCC,LOW,HIGH
      DIMENSION EMAT(NPARM,NPARM)
      CHARACTER*10 CHNAM
      DOUBLE PRECISION VAL,ERROR,BND1,BND2
      INTEGER IVARBL
      INTEGER NSTEPS 
      DOUBLE PRECISION LEFT,RIGHT,VAR,INTEGRAL,STEP
      IF (IFLAG.NE.1) GOTO 100
      
 100  SUM=1.D0
      NSTEPS = 1000
      LEFT  = 1.5318D0-5.6D-3*10.D0
      RIGHT = 1.5318D0+5.6D-3*10.D0
      STEP  = (RIGHT-LEFT)/DFLOAT(NSTEPS)
      INTEGRAL = 0.D0
      DO I=1,NSTEPS
        VAR = LEFT + (DFLOAT(I)-0.5D0)*STEP
        INTEGRAL = INTEGRAL + BWG(VAR,U)*STEP
      ENDDO
      SUM = INTEGRAL
      

C     
C     Gaussian Constraints
C     
      GCG   = (U(3)-9.1D-3)/0.5D-3
      GCG = 0.D0
      
      DO I=1,NDATA
        IF (DABS((XDATA(I)-1.5318D0)/5.6D-3).LE.10.D0) THEN
          FUN= BWG(XDATA(I),U)
          IF (FUN.LE.0) THEN
            SUM = SUM + 1000.0
          ELSE
            SUM=SUM-DLOG(FUN)
          ENDIF
        ENDIF
      END DO
      SUM = SUM + 0.5*GCG*GCG
      F=SUM*2.D0
      IF (IFLAG.EQ.3) THEN
C     WRITE(22,1000  )F
C     1000   FORMAT( G20.10 )
        DO I=1,NPARM
C     CALL MNERRS(-I,EPLUS,EMINUS,EPARAB,GLOBCC)
          CALL MNPOUT (I,CHNAM,VAL,ERROR,BND1,BND2,IVARBL)
          EU(I)=ERROR
          WRITE(22,1003  )I,CHNAM,VAL,ERROR,BND1,BND2,IVARBL
 1003     FORMAT(7X,I5,1X,A10,4G20.10,I5)
        ENDDO
        WRITE(23,1004)(U(3)-9.1D-3)/EU(3),(U(4)-1.5318D0)/EU(4),
     &       (U(1)-2000.D0)/EU(1)
 1004   FORMAT (3G20.10)
      ENDIF
      RETURN
      END
      
      
 
      FUNCTION GAUSS(OBM,TRM,SIG)
      IMPLICIT NONE 
      DOUBLE PRECISION X,OBM,TRM,SIG,GAUSS,DSQRT,DEXP
      X = (OBM - TRM) / SIG
      GAUSS = 1.D0/DSQRT(2.D0*3.14159265358979312D0)/SIG*
     &     DEXP(-0.5D0*X**2)
      RETURN
      END
C     
C     homemade function
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

      NPOINTS=100
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
      EXTERNAL WWERF
      COMPLEX*16 Z,WZ,WWERF
      DOUBLE PRECISION X,Y,FOLD,OBM,GAM,SIG,TRM,ARG
      ARG = (OBM-1.5318D0)/5.6D-3
      IF (DABS(ARG)>10.D0) THEN 
        FOLD = 0.D0
        RETURN
      ENDIF
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
      DOUBLE PRECISION N2,SIG2,M2,RATIO,SIG3,PWR
      DOUBLE PRECISION BWG,FOLD,DSQRT,GAUSS,FOLD1
      DIMENSION U(*)
      
      
      
      N1       = U(1)
      SIG1     = U(2)
      GAMMA    = U(3)
      M1       = U(4)
C
C 9.48388017490093271D-01   = int_10^10 BW dm  
C
      
      BWG=N1*FOLD(X,GAMMA,SIG1,M1)
      
      RETURN
      END

      FUNCTION BWG2(X,EX,U)
      IMPLICIT NONE
      DOUBLE PRECISION X,EX,U
      DOUBLE PRECISION N1,SIG1,GAMMA,M1,A0,A1,A2,A3,DM
      DOUBLE PRECISION N2,SIG2,M2,RATIO,SIG3,PWR
      DOUBLE PRECISION BWG2,FOLD,DSQRT,GAUSS,FOLD1
      DIMENSION U(*)
      


      N1       = U(1)
      SIG1     = U(2)*EX
      GAMMA    = U(3)
      M1       = U(4)

      PWR      = U(5)
      A0       = U(6)
      A1       = U(7)
      A2       = U(8)
      A3       = U(9)

      DM = X - (1.321D0+0.1396D0)

      IF (DM.LT.0) THEN
        BWG2 = 0.D0
        RETURN
      ENDIF

C
C 9.36548965138892964D-01 = int_10^10 BW dm  
C
      BWG2=N1*FOLD(X,GAMMA,SIG1,M1)/9.36548965138892964D-01+
     &     (DM**PWR)*
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

      FUNCTION WIDTH1(X)
      DOUBLE PRECISION X
      DOUBLE PRECISION WIDTH1
      WIDTH1 = 9.17220D-5*X+1.20072D-3 
      RETURN
      END

      FUNCTION WIDTH2(X)
      DOUBLE PRECISION X
      DOUBLE PRECISION WIDTH1
      WIDTH1 = 2.06315D-04*X+2.82050D-3
      RETURN
      END


  
      SUBROUTINE PRTERR 
C   a little hand-made routine to print out parameter errors 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
C  find out how many variable parameters there are 
      CALL MNSTAT(FMIN,FEDM,ERRDEF,NPARI,NPARX,ISTAT) 
C   and their errors 
      DO 50 I= 1, NPARI 
      CALL MNERRS(-I,EPLUS,EMINUS,EPARAB,GLOBCC) 
      WRITE (6,45) I,EPLUS,EMINUS,EPARAB,GLOBCC 
   45 FORMAT (5X,I5,4F12.6) 
   50 CONTINUE 
      RETURN 
      END 
