C*
C* this one to fit xi- pi- with just background function
C*
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
      REAL RNDM,DUMMY
      EXTERNAL FCN
      PARAMETER(NPARM=5)
      CHARACTER*10 PNAM
      INTEGER NPRM
      DOUBLE PRECISION VSTRT, STP, ARGLIS
      COMMON/PARS/ VSTRT(NPARM), STP(NPARM),ARGLIS(10),
     &     NPRM(NPARM), PNAM(NPARM) 

      INTEGER IERFLG
      DATA NPRM /   1,2,3,4,5/
      DATA PNAM / 'PWR','A0','A1','A2','A3'/
      DATA VSTRT/ 
     &             0.55732  , 0.32449D+06, -0.80838D+06, 
     &             0.10004D+07, -0.45648D+06/
      DATA STP  /   5*0.001 /
C      DATA STP  /   5*0.0001, 0.000D0, 0.D0, 0.D0 /


      CALL HLIMIT(300000)
      CALL HROPEN(20,'blah','signal.hist','N',1024,ISTAT)
      CALL HBOOK1(1,'mass',280,1.4,2.8,0.0)


      NDATA=1
C      OPEN(UNIT=80,FILE='data_points_wrong_sign.dat')
      OPEN(UNIT=80,FILE='jet20_data_points_ws.dat')
C      OPEN(UNIT=80,FILE='data_points_unique_ws.dat')
C      OPEN(UNIT=80,FILE='data_points_ws_deltaR.dat')
 50   READ(80,*,END=500)NUMBER,RUN,EVENT,
     &     XDATA(NDATA),EXDATA(NDATA)
      CALL HF1(1,REAL(XDATA(NDATA)),1.0)
      NDATA=NDATA+1
      GOTO 50
 500  NDATA=NDATA-1
      GOTO 999
 800  FORMAT(1X,'Total data pts read:  ',i10)
 999  WRITE(*,800)NDATA
      CLOSE(80)

      CALL HROUT (0,ICYCLE,' ')
      CALL HREND('blah')

C*
C*   Open data files
C*
      open(unit=21,file='fcn_xipiminus_only1.dat',
     &     status='unknown',form='formatted')
      open(unit=22,file='func_xipiminus_only1.dat',
     &     status='unknown',form='formatted')
      open(unit=23,file='emat_xipiminus_only1.dat',
     &     status='unknown',form='formatted')

      CALL MNINIT(5,6,7) 

      DO I=1,NPARM
        CALL MNPARM(NPRM(I),PNAM(I),VSTRT(I),STP(I),0.D0,0.D0,IERFLG) 
        IF (IERFLG .NE. 0)  THEN 
          WRITE (6,'(A,I)')  ' UNABLE TO DEFINE PARAMETER NO.',I 
          STOP 
        ENDIF  
      ENDDO


      CALL MNSETI('xicst')

      
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


      CLOSE(5)
      CLOSE(21)
      CLOSE(22)

      STOP
      END


      SUBROUTINE FCN(NPAR,GIN,F,U,IFLAG)
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER IFLAG,NPAR
      DOUBLE PRECISION U,GIN,F
      DIMENSION U(*),GIN(*)
      INTEGER MXDATA,NDATA,NPARM,I,NPOINTS
      PARAMETER(MXDATA=1000000)
      DOUBLE PRECISION XDATA,EXDATA
      COMMON /USER/ XDATA(MXDATA),EXDATA(MXDATA),NDATA
      DOUBLE PRECISION FUN,SUM,SQRT,DEXP,BWG,BWG2,C1,C2,HIGHM
      DOUBLE PRECISION GCM,GCS,GCG
      PARAMETER(NPARM=5)
      CHARACTER*10 PNAM
      INTEGER NPRM
      DOUBLE PRECISION VSTRT, STP, ARGLIS
      COMMON/PARS/ VSTRT(NPARM), STP(NPARM),ARGLIS(10),
     &     NPRM(NPARM), PNAM(NPARM)
      DOUBLE PRECISION EU,FOLD
      DIMENSION EU(NPARM)
      DOUBLE PRECISION EMAT,EPLUS,EMINUS,EPARAB,GLOBCC,LOW,HIGH
      DIMENSION EMAT(NPARM,NPARM)
      CHARACTER*10 CHNAM
      DOUBLE PRECISION VAL,ERROR,BND1,BND2
      INTEGER IVARBL
      INTEGER NSTEPS,J
      DOUBLE PRECISION LEFT,RIGHT,VAR,STEP
      INTEGER NPARI,NPARX,ISTAT
      DOUBLE PRECISION FMIN,FEDM,ERRDEF 
      DOUBLE PRECISION SIG,DSIG,MASS



      HIGHM = 2.4D0
      C1    = HIGHM - (1.321D0+0.1396D0)

      IF (IFLAG.NE.1) GOTO 100
      
 100  SUM=0.D0

      SUM =C1**(U(1))*(U(2)*C1/(1.D0+U(1))+
     &     U(3)*C1**2/(2.D0+U(1))+
     &     U(4)*C1**3/(3.D0+U(1))+
     &     U(5)*C1**4/(4.D0+U(1)))     


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
      SUM = SUM
      F=SUM*2.D0
      IF (IFLAG.EQ.3) THEN
        WRITE(21,1000  )F
 1000   FORMAT( G20.10 )
        DO I=1,NPARM
          CALL MNPOUT (I,CHNAM,VAL,ERROR,BND1,BND2,IVARBL)
          WRITE(22,1003)I,CHNAM,VAL,ERROR,BND1,BND2,IVARBL
 1003     FORMAT(7X,I5,1X,A10,4G20.10,I5)
        ENDDO
C*
C*  find out how many variable parameters there are
C*
        
          CALL MNSTAT(FMIN,FEDM,ERRDEF,NPARI,NPARX,ISTAT)
C*
C* write out error matrix
C*
        CALL MNEMAT(EMAT,NPARM)
        DO I=1,NPARI
           WRITE(23,'(2X,13G20.10)')(EMAT(I,J),J=1,NPARI)
        ENDDO

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
      DOUBLE PRECISION ARG
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
      DOUBLE PRECISION LEFT,RIGHT,VAR,INTEGRAL,STEP,
     &     CONST
      INTEGER NSTEPS,I
      

      PWR      = U(1)
      A0       = U(2)
      A1       = U(3)
      A2       = U(4)
      A3       = U(5)
C      N1       = U(6)
C      SIG1     = U(7)
C      M1       = U(8)
C
      DM = X - (1.321D0+0.1396D0)

      IF (DM.LT.0) THEN
        BWG = 0.D0
        RETURN
      ENDIF

      BWG=
     &     (DM**PWR)*
     &     (A0+A1*DM+A2*DM**2+A3*DM**3)
      
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

      BWG2=N1*FOLD(X,GAMMA,SIG1,M1)+
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

C
C this subroutine calculates width
C      
      
      SUBROUTINE WIDTH(X,SIG,DSIG)
      DOUBLE PRECISION X, MASS,SIG,DSIG
      DOUBLE PRECISION EMAT, PAR, EPAR
      DIMENSION EMAT(2,2), PAR(2), EPAR(2)
      
      DATA PAR / 0.9380063955D-02, -10.85937791D0 / 
      DATA EPAR / 0.2115535481D-02, 3.447891888D0 /

      DATA EMAT / 0.4475490371D-05, -0.7201371655D-02, 
     &     -0.7201371655D-02, 11.88795847D0 /

      MASS = 1000.D0*X
      SIG = PAR(1)*MASS + PAR(2)
      SIG = 1.D-3 * SIG
      DSIG = DSQRT(MASS*MASS*EMAT(1,1)+EMAT(2,2)
     &     + 2.D0*MASS*EMAT(1,2))
      DSIG = 1.D-3 * DSIG
                  
      RETURN 
      END

C
C this subroutine calculates efficiency
C

      SUBROUTINE EFFIC(X,E,DE)
      DOUBLE PRECISION X, E,DE,MASS
      DOUBLE PRECISION  PAR, EPAR
      
      DATA PAR  / 0.6506024096D-03 /
      DATA EPAR / 0.1077623122E-03 / 

      MASS = 1.D+3 * X

      E  = PAR * ( MASS - 1531.D0) + 1.D0 
      DE = EPAR * ( MASS  - 1531.D0) 
                  
      RETURN 
      END


