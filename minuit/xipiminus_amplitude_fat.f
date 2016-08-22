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

      NDATA=1
      OPEN(UNIT=80,FILE='data_points_wrong_sign.dat')
C      OPEN(UNIT=80,FILE='jet20_data_points_ws.dat')
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
      open(unit=5,file='xipiminus_amplitude.data',status='old')
      open(unit=21,file='fcn_xipiminus_amplitude_fat.dat',
     &     status='unknown',form='formatted')
      open(unit=22,file='func_xipiminus_amplitude_fat.dat',
     &     status='unknown',form='formatted')
      open(unit=23,file='emat_xipiminus_amplitude_fat.dat',
     &     status='unknown',form='formatted')
      open(unit=24,file='ul_xipiminus_amplitude_fat.dat',
     &     status='unknown',form='formatted')
C      open(unit=6,file='fit.out',status='unknown',form='formatted')

      CALL MINTIO(5,6,7)
      CALL MINUIT(FCN,0)

      CLOSE(5)
      CLOSE(21)
      CLOSE(22)
      CLOSE(23)
      CLOSE(24)

      STOP
      END


      SUBROUTINE FCN(NPAR,GIN,F,U,IFLAG)
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
      DOUBLE PRECISION GCM,GCS,GCG,GCE,CONSTR,GCS1,GCG1
      INTEGER NPARM
      PARAMETER(NPARM=9)
      DOUBLE PRECISION EMAT,EPLUS,EMINUS,EPARAB,GLOBCC,MASS
      DIMENSION EMAT(NPARM,NPARM)
      CHARACTER*10 CHNAM
      INTEGER IVARBL
      DOUBLE PRECISION VAL,ERROR,BND1,BND2
      INTEGER I,J,I90,I95,HAS90,HAS95
      DOUBLE PRECISION R90, R9 
      DIMENSION LH(1001)
      DOUBLE PRECISION LHSUM,LHRSUM,LH,RINT,RF
      DOUBLE PRECISION EFF,DEFF,SIG,DSIG
      DOUBLE PRECISION LEFT,RIGHT,VAR,INTEGRAL,STEP,INTEGRAL1
      DOUBLE PRECISION ARG,SIGMA,DSIGMA,SIGAM
      INTEGER NSTEPS
      DOUBLE PRECISION FOLD,FOLDER
      DOUBLE PRECISION FMIN,FEDM,ERRDEF 
      INTEGER NPARI,NPARX,ISTAT


      
      CALL VZERO(LH,1001)

      HIGHM = 2.2D0
      C1    = HIGHM - (1.321D0+0.1396D0)

      
      MASS  = 1.862D0
      CALL WIDTH(MASS,SIG,DSIG)
      CALL EFFIC(MASS,EFF,DEFF)


C*
C*  introduce systematic error on sigmas (0.5 MeV)  and efficiencies (0.1)
C* 

      DEFF = DSQRT(DEFF*DEFF+1.D-2*EFF*EFF)
      DSIG = DSQRT(DSIG*DSIG+0.5D-3*0.5D-3)

      IF (IFLAG.NE.1) GOTO 100
      
 100  SUM=0.D0
 
      NSTEPS = 100


      CALL WIDTH(1.862D0,SIGMA,DSIGMA)
C*
C*  introduce systematic error on sigmas (0.5 MeV)  and efficiencies (0.1)
C* 
      DSIGMA = DSQRT(DSIGMA*DSIGMA+0.5D-3*0.5D-3)
      SIGAM = DSQRT(SIGMA*SIGMA+0.25D0*0.017D0*0.017D0)

      LEFT  = 1.862D0-SIGAM*10.D0
      RIGHT = 1.862D0+SIGAM*10.D0
      STEP  = (RIGHT-LEFT)/DFLOAT(NSTEPS)

      INTEGRAL = 0.D0
      DO I=1,NSTEPS
        VAR = LEFT + (DFLOAT(I)-0.5D0)*STEP
        INTEGRAL = INTEGRAL+FOLD(VAR,U(3),U(2),U(4))*STEP
      ENDDO
      
      SUM=U(1)*INTEGRAL+C1**(U(5))*(U(6)*C1/(1.D0+U(5))+
     &     U(7)*C1**2/(2.D0+U(5))+
     &     U(8)*C1**3/(3.D0+U(5))+
     &     U(9)*C1**4/(4.D0+U(5)))
C
C     Gaussian Constraints
C     
      GCG1   = (U(3)-3.6D-3)/2.0D-3
      GCS    = (U(2)-SIG)/DSIG
      GCM    = (U(4)-MASS)/2.0D-3
 
      CONSTR = 0.5*GCS*GCS+0.5*GCM*GCM+0.5*GCE*GCE+
     &     0.5*GCG1*GCG1
 

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
      SUM = SUM + CONSTR
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
C*
C* doing LH scan 
C*
        PRINT*,'Starting LH scan'
        RINT   = 0.D0
        LHSUM  = 0.D0
        LHRSUM = 0.D0
        HAS90  = 0
        HAS95  = 0
        I90    = 0
        I95    = 0

        DO I=1,1001
          U(1)  = DFLOAT(I-1)
          SUM=U(1)*INTEGRAL+C1**(U(5))*(U(6)*C1/(1.D0+U(5))+
     &         U(7)*C1**2/(2.D0+U(5))+
     &         U(8)*C1**3/(3.D0+U(5))+
     &         U(9)*C1**4/(4.D0+U(5)))
          
          DO IDATA=1,NDATA
            IF (XDATA(IDATA).LE.HIGHM) THEN
              FUN= BWG(XDATA(IDATA),U)
              IF (FUN.LE.0) THEN
                SUM = SUM + 1000.0
              ELSE
                SUM=SUM-DLOG(FUN)
              ENDIF
            ENDIF
          END DO
          SUM = SUM + CONSTR
          RF=SUM*2.D0
          LH(I)  = DEXP(-0.5D0*(RF-F))
          LHSUM  = LHSUM + LH(I)
        ENDDO
        DO I=1,1001

          WRITE(21,1004  )LH(I)
 1004     FORMAT( G20.10 )

          LHRSUM = LHRSUM + LH(I)
          IF ( LHRSUM / LHSUM .GE. 0.9. AND. HAS90.EQ.0 ) THEN
            I90 = I
            HAS90=1
          ENDIF
          IF ( LHRSUM / LHSUM .GE. 0.95. AND. HAS95.EQ.0 ) THEN
            I95 = I
            HAS95=1
          ENDIF
        ENDDO

       WRITE(24,1003  )I90,I95,DFLOAT(I90-1)*INTEGRAL,
     &      DFLOAT(I95-1)*INTEGRAL
       
 1003  FORMAT(2X,2I5)
        
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
      DOUBLE PRECISION ARG,SIGMA,DSIGMA,SIGAM


      CALL WIDTH(1.862D0,SIGMA,DSIGMA)
      SIGAM = DSQRT(SIGMA*SIGMA+0.25D0*0.017D0*0.017D0)
      ARG = (OBM-1.862D0)/SIGAM

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
      DOUBLE PRECISION N2,SIG2,M2,RATIO,SIG3,PWR,GAM2
      DOUBLE PRECISION BWG,FOLD,DSQRT,GAUSS,FOLD1
      DIMENSION U(*)
      


      N1       = U(1)
      SIG1     = U(2)
      GAMMA    = U(3)
      M1       = U(4)

      PWR      = U(5)
      A0       = U(6)
      A1       = U(7)
      A2       = U(8)
      A3       = U(9)


      DM = X - (1.321D0+0.1396D0)

      IF (DM.LT.0) THEN
        BWG = 0.D0
        RETURN
      ENDIF

      BWG=N1*FOLD(X,GAMMA,SIG1,M1)+
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

