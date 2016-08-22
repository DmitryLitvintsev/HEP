C
C this is toy MC generation to test coverage, Xi pi-
C 
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
      REAL RNDM2,DUMMY,RNDM,RAND
      EXTERNAL FCN
      DOUBLE PRECISION BWG

      PARAMETER(NPARM=9)
      CHARACTER*10 PNAM
      INTEGER NPRM
      DOUBLE PRECISION VSTRT, STP, ARGLIS
      COMMON/PARS/ VSTRT(NPARM), STP(NPARM),ARGLIS(10),
     &     NPRM(NPARM), PNAM(NPARM) 



      INTEGER IERFLG
      DATA NPRM /   1,2,3,4,5,6,7,8,9  / 
      DATA PNAM / 'N1', 'SIG1', 'GAMMA','M1','PWR',
     &            'A0','A1','A2','A3' /

      DATA VSTRT/  100.D0  , 6.5D-3, 0.D0, 1.862D0, 
     &             0.6481679573D0, 332112.4879D0, 
     &             -836553.1246D0, 994755.1096D0,
     &             -422324.6490D0 /
      DATA STP  /   0.001 ,    0.0001 ,  0.00,   6*0.001 /



      CALL HLIMIT(300000)
      CALL HROPEN(20,'blah','signal1.hist','N',1024,ISTAT)
      CALL HBOOK1(1,'mass',280,1.4,2.8,0.0)

      CALL SRAND(time())

      PRINT *, RAND()


      MASS0  =  1.862D0
C      NDATA  = 103537 
C    dtaa points with mass < 2.2 

      NDATA  = 48815 
      NEXP   = 10
      SIGMA  = 3.2D-3
      GAMMA  = 9.1D-3
      IDATA  = 0
      I      = 0

      
C*     
C*   Open data files
C*     
C      open(unit=5,file='signal_only.data',status='old')
      open(unit=22,file='parms_minus.dat',
     &     status='unknown',form='formatted')
      open(unit=23,file='residuals_minus.dat',
     &     status='unknown',form='formatted')
      open(unit=24,file='ul_xipipinus_coverage.dat',
     &     status='unknown',form='formatted')
C      open(unit=6,file='parms_minus.out',
C     &     status='unknown',form='formatted')
      


      J = 0 
      

      CALL MNINIT(5,6,7) 
      CALL MNSETI('xicst')
      
 1    CONTINUE


      DO I=1,NPARM
        CALL MNPARM(NPRM(I),PNAM(I),VSTRT(I),STP(I),0.D0,0.D0,IERFLG) 
        IF (IERFLG .NE. 0)  THEN 
          WRITE (6,'(A,I)')  ' UNABLE TO DEFINE PARAMETER NO.',I 
          STOP 
        ENDIF  
      ENDDO

      IDATA = 0
      DO WHILE ( IDATA .LT. NDATA ) 
        MASS = 1.4D0 + RAND()*(2.8D0-1.4D0)
        IF (MASS.LT.2.2D0) THEN 
          PRB = 1000000.D0*RAND()
          IF (PRB.LE.BWG(MASS,VSTRT)) THEN
            IDATA = IDATA+1
            XDATA(IDATA) = MASS 
            CALL HF1(1,REAL(MASS),1.0)
          ENDIF
        ENDIF
      ENDDO       
      
        
      ARGLIS(1) = 1. 
      CALL MNEXCM(FCN, 'CALL FCN', ARGLIS ,1,IERFLG,0) 

      ARGLIS(1) = 2. 
C     CALL MNEXCM(FCN,'FIX',       ARGLIS ,1,IERFLG,0) 
      CALL MNEXCM(FCN,'SET PRINT', ARGLIS ,1,IERFLG,0) 
C      CALL MNEXCM(FCN,'SET STR',   ARGLIS ,1,IERFLG,0) 
      ARGLIS(1) = 0.
      CALL MNEXCM(FCN,'MIGRAD', ARGLIS ,0,IERFLG,0) 
C      CALL MNEXCM(FCN,'MINOS', ARGLIS ,0,IERFLG,0) 
      ARGLIS(1) = 3.
      CALL MNEXCM(FCN,'CALL FCN', ARGLIS , 1,IERFLG,0) 
      CALL MNEXCM(FCN,'STOP ', 0.D0,0,IERFLG,0)

C      CALL MINUIT(FCN,0)
      J=J+1
      IF (J.LT.NEXP) GOTO 1

      CALL HROUT (0,ICYCLE,' ')
      CALL HREND('blah')

      
      CLOSE(5)
C     CLOSE(6)
      CLOSE(22)
      CLOSE(23)
      CLOSE(24)
      
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
      DOUBLE PRECISION GCM,GCS,GCG,CONSTR
      PARAMETER(NPARM=9)
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
      DOUBLE PRECISION LEFT,RIGHT,VAR,INTEGRAL,STEP
      DOUBLE PRECISION SIG,DSIG,MASS

      INTEGER I90,I95,HAS90,HAS95,IDATA
      DOUBLE PRECISION R90, R9 
      DIMENSION LH(1001)
      DOUBLE PRECISION LHSUM,LHRSUM,LH,RINT,RF


      
      HIGHM = 2.2D0
      C1    = HIGHM - (1.321D0+0.1396D0)
      
      MASS  = 1.862D0
      CALL WIDTH(MASS,SIG,DSIG)

      
      IF (IFLAG.NE.1) GOTO 100
      
 100  SUM=1.D0
      INTEGRAL = 1.D0
      SUM = U(1)+
     &     C1**(U(5))*(U(6)*C1/(1.D0+U(5))+
     &     U(7)*C1**2/(2.D0+U(5))+
     &     U(8)*C1**3/(3.D0+U(5))+
     &     U(9)*C1**4/(4.D0+U(5)))
C     
C     Gaussian Constraints
C     
      GCM   = (U(4)-1.862D0)/2.0D-3
      GCS   = (U(2)-SIG)/DSIG

      CONSTR = 0.5*GCS*GCS+0.5*GCM*GCM

      DO I=1,NDATA
        IF (XDATA(I).LT.HIGHM) THEN
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
C     WRITE(22,1000  )F
C     1000   FORMAT( G20.10 )
        DO I=1,NPARM
C     CALL MNERRS(-I,EPLUS,EMINUS,EPARAB,GLOBCC)
          CALL MNPOUT (I,CHNAM,VAL,ERROR,BND1,BND2,IVARBL)
          EU(I)=ERROR
          IF (I.EQ.1) THEN
            VAL   = VAL * INTEGRAL
            ERROR = ERROR * INTEGRAL 
            U(I)  = U(I) * INTEGRAL
            EU(I) = EU(I) * INTEGRAL
          ENDIF
          IF (I.EQ.1) THEN
            WRITE(22,1003  )I,CHNAM,VAL,ERROR,BND1,BND2,IVARBL
 1003       FORMAT(7X,I5,1X,A10,4G20.10,I5)
          ENDIF
        ENDDO
        WRITE(23,1004)((U(J)-VSTRT(J))/EU(J),J=1,NPARM)
 1004   FORMAT (9G20.10)


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
          U(1)  =  DFLOAT(I-1)
          SUM=U(1)+C1**(U(5))*(U(6)*C1/(1.D0+U(5))+
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
          SUM = SUM + 0.5*GCS*GCS + 0.5*GCM*GCM
          RF=SUM*2.D0
          LH(I)  = DEXP(-0.5D0*(RF-F))
          LHSUM  = LHSUM + LH(I)
        ENDDO

        DO I=1,1001
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
        
       WRITE(24,1001  )I90-1,I95-1
 1001  FORMAT(2X,2I5)

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

      IF (DM.LT.0D0) THEN
        BWG = 0.D0
        RETURN
      ENDIF

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
