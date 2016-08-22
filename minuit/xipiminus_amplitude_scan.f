C
C unbinned likelihood fit  of Xi- pi-
C
      PROGRAM MAIN
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT NONE 
      INTEGER MXDATA,NDATA
      PARAMETER(MXDATA=1000000)
      DOUBLE PRECISION XDATA,EXDATA
      COMMON /USER/ XDATA(MXDATA),EXDATA(MXDATA),NDATA
      INTEGER ISTAT,NUMBER,RUN,EVENT
      REAL RNDM,DUMMY
      EXTERNAL FCN

      INTEGER NPARM,NEXP,I,J,K
      PARAMETER(NPARM=8)
      CHARACTER*10 PNAM
      INTEGER NPRM
      DOUBLE PRECISION VSTRT, STP, ARGLIS
      COMMON/PARS/ VSTRT(NPARM), STP(NPARM),ARGLIS(10),
     &     NPRM(NPARM), PNAM(NPARM) 



      INTEGER IERFLG
      DATA NPRM /   1,2,3,4,5,6,7,8  / 
      DATA PNAM / 'N1', 'SIG1','M1','PWR',
     &            'A0','A1','A2','A3' /

C     ttt
      DATA VSTRT/  10.D0  , 6.5D-3, 1.862D0, 
     &              0.64817D0, 0.33211D+06, 
     &             -0.83655D+06, 0.99476D+06, -0.42232D+06 /
C jet20

C      DATA VSTRT/  0.D0  , 6.5D-3, 1.862D0, 
C     &             0.58214D0, 41597.D0,
C     &             -0.10206D+06, 0.11786D+06, -49349.D0/

      DATA STP  /   0.001D0 , 0.1D-3, 0.01D0,  0.01D0,  4*0.001 / 

      
      OPEN(UNIT=10,FILE='number.txt')
      READ(10,*,END=100)NEXP
 100  CONTINUE
      PRINT*,NEXP 


      NDATA=1
C      OPEN(UNIT=80,FILE='data_points_wrong_sign.dat')
      OPEN(UNIT=80,FILE='jet20_data_points_ws.dat')
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
C*   Open data files
C*
      open(unit=21,file='fcn_xipiminus_amplitude_scan.dat',
     &     status='unknown',form='formatted')
      open(unit=22,file='func_xipiminus_amplitude_scan.dat',
     &     status='unknown',form='formatted')
      open(unit=23,file='emat_xipiminus_amplitude_scan.dat',
     &     status='unknown',form='formatted')
      open(unit=24,file='ul_xipiminus_amplitude_scan.dat',
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
      CALL MNEXCM(FCN,'SET PRINT', ARGLIS ,1,IERFLG,0) 
      ARGLIS(1) = 2. 
      CALL MNEXCM(FCN,'SET STR', ARGLIS ,1,IERFLG,0) 
      ARGLIS(1) = 2. 
      CALL MNEXCM(FCN,'FIX PAR', ARGLIS ,1,IERFLG,0) 
      ARGLIS(1) = 3. 
      CALL MNEXCM(FCN,'FIX PAR', ARGLIS ,1,IERFLG,0) 
      CALL MNEXCM(FCN,'MIGRAD', ARGLIS ,0,IERFLG,0) 
      ARGLIS(1) = 2. 
      CALL MNEXCM(FCN,'REL PAR', ARGLIS ,1,IERFLG,0) 
      ARGLIS(1) = 3. 
      CALL MNEXCM(FCN,'REL PAR', ARGLIS ,1,IERFLG,0) 
      

C      I = (NEXP-1)*10
 1    CONTINUE

      ARGLIS(1) = 1.
      ARGLIS(2) = DFLOAT(I)
      CALL MNEXCM(FCN,'SET PAR ', ARGLIS ,2,IERFLG,0) 
      CALL MNEXCM(FCN,'MIGRAD', ARGLIS ,0,IERFLG,0) 
      ARGLIS(1) = 3.
      CALL MNEXCM(FCN,'CALL FCN', ARGLIS , 1,IERFLG,0) 
      I=I+1
C      IF (I.LT.NEXP*10) GOTO 1
      CALL MNEXCM(FCN,'STOP ', 0.D0,0,IERFLG,0)

      CLOSE(21)
      CLOSE(22)
      CLOSE(23)
      STOP
      END


      SUBROUTINE FCN(NPAR,GIN,F,U,IFLAG)
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT NONE 
      DOUBLE PRECISION GIN,F,U
      DIMENSION U(*),GIN(*)
      INTEGER NPAR,IFLAG
      INTEGER MXDATA,NDATA
      PARAMETER(MXDATA=1000000)
      DOUBLE PRECISION XDATA,EXDATA
      COMMON /USER/ XDATA(MXDATA),EXDATA(MXDATA),NDATA
      INTEGER IDATA
      DOUBLE PRECISION FUN,SUM,SQRT,DEXP,BWG,BWG2,C1,C2,HIGHM
      DOUBLE PRECISION GCM,GCS,GCG
      INTEGER NPARM
      PARAMETER(NPARM=9)
      DOUBLE PRECISION EMAT,EPLUS,EMINUS,EPARAB,GLOBCC
      DIMENSION EMAT(NPARM,NPARM)
      INTEGER I,J,I90,I95,HAS90,HAS95
      DOUBLE PRECISION R90, R9 
      DIMENSION LH(1001)
      DOUBLE PRECISION LHSUM,LHRSUM,LH,RINT,RF

      CHARACTER*10 CHNAM
      INTEGER IVARBL
      DOUBLE PRECISION VAL,ERROR,BND1,BND2,MASS

      DOUBLE PRECISION EFF,DEFF,SIG,DSIG


      HIGHM = 2.4D0
      C1    = HIGHM - (1.321D0+0.1396D0)

      MASS  = 1.862D0
      CALL WIDTH(MASS,SIG,DSIG)
      CALL EFFIC(MASS,EFF,DEFF)
      
      IF (IFLAG.NE.1) GOTO 100
      
 100  SUM=U(1)+C1**(U(4))*(U(5)*C1/(1.D0+U(4))+
     &     U(6)*C1**2/(2.D0+U(4))+
     &     U(7)*C1**3/(3.D0+U(4))+
     &     U(8)*C1**4/(4.D0+U(4)))

C
C Gaussian Constraints
C
      GCS    = (U(2)-SIG)/DSIG
      GCM    = (U(3)-1.862D0)/2.0D-3

C      GCS   = 0.D0

      DO I=1,NDATA
        IF (XDATA(I).LE.HIGHM) THEN
          FUN= BWG(XDATA(I),U)
C          FUN= BWG2(XDATA(I),EXDATA(I),U)
          IF (FUN.LE.0) THEN
            SUM = SUM + 1000.0
          ELSE
            SUM=SUM-DLOG(FUN)
          ENDIF
        ENDIF
      END DO
      SUM = SUM + 0.5D0*GCS*GCS + 0.5D0*GCM*GCM
C      SUM = SUM  +  + 0.5D0*GCM*GCM
      F=SUM*2.D0
      IF (IFLAG.EQ.3) THEN
C*
C* write out FCN
C*
        WRITE(21,1000  )F
 1000   FORMAT( G20.10 )

C*
C* write out function parameters
C*
        DO I=1,NPARM
          CALL MNERRS(-I,EPLUS,EMINUS,EPARAB,GLOBCC)
          WRITE(22,1001  )I,U(I),EPARAB,EPLUS,EMINUS,GLOBCC
 1001     FORMAT(6X,I5,5G20.10)
        ENDDO
C*
C* write out error matrix
C*
        CALL MNEMAT(EMAT,NPARM)
        DO I=1,NPARM
          WRITE(23,1002)(EMAT(I,J),J=1,NPARM)
 1002     FORMAT(8x,8G20.10)
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
      DOUBLE PRECISION N2,SIG2,M2,RATIO,SIG3,PWR,GAM2
      DOUBLE PRECISION BWG,FOLD,DSQRT,GAUSS,FOLD1
      DIMENSION U(*)
      


      N1       = U(1)
      SIG1     = U(2)
      M1       = U(3)

      PWR      = U(4)
      A0       = U(5)
      A1       = U(6)
      A2       = U(7)
      A3       = U(8)
      


      DM = X - (1.321D0+0.1396D0)

      IF (DM.LT.0) THEN
        BWG = 0.D0
        RETURN
      ENDIF

C      BWG=N1*FOLD(X,GAMMA,SIG1,M1)+

      BWG=N1*GAUSS(X,M1,SIG1)+
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

