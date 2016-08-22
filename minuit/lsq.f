C
C unbinned likelihood fit  of Xi(1530)  only
C
      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER MXDATA,NDATA,NPARM
      PARAMETER(MXDATA=1000000)
      DOUBLE PRECISION X,Y,DY
      COMMON /USER/ X(MXDATA),Y(MXDATA),DY(MXDATA),NDATA
      INTEGER I, J, IDATA, ICYCLE
      EXTERNAL FCN

      PARAMETER(NPARM=2)
      CHARACTER*10 PNAM
      INTEGER NPRM
      DOUBLE PRECISION VSTRT, STP, ARGLIS
      COMMON/PARS/ VSTRT(NPARM), STP(NPARM),ARGLIS(10),
     &     NPRM(NPARM), PNAM(NPARM) 

      INTEGER IERFLG
      DATA NPRM /   1 , 2 /
      DATA PNAM / 'K', 'B' /
      DATA VSTRT/ 9.D-3, -11.D0/
      DATA STP  /   0.0001 ,    0.0001/

      NDATA=2
      X(1) = 1531.D0
      X(2) = 2469.D0
      Y(1) = 3.5015D0
      DY(1) = 0.57248D0
      Y(2) = 12.3D0
      DY(2) = 1.9D0

      open(unit=22,file='lsq.dat',
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
      CALL MNEXCM(FCN,'SET PRINT', ARGLIS ,1,IERFLG,0) 
      CALL MNEXCM(FCN,'SET STR', ARGLIS ,1,IERFLG,0) 
      CALL MNEXCM(FCN,'MIGRAD', ARGLIS ,0,IERFLG,0) 
      ARGLIS(1) = 3.
      CALL MNEXCM(FCN,'CALL FCN', ARGLIS , 1,IERFLG,0) 
      CALL MNEXCM(FCN,'STOP ', 0.D0,0,IERFLG,0)

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
      DOUBLE PRECISION X,Y,DY
      COMMON /USER/ X(MXDATA),Y(MXDATA),DY(MXDATA),NDATA
      DOUBLE PRECISION FUN,SUM,SQRT,DEXP,BWG,BWG2,C1,C2,HIGHM
      DOUBLE PRECISION GCM,GCS,GCG
      PARAMETER(NPARM=2)
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


 100  SUM=0.D0
      DO I=1,NDATA
        SUM = SUM + (U(1)*X(I)+U(2)-Y(I))**2/DY(I)/DY(I)
      END DO
      F=SUM
      IF (IFLAG.EQ.3) THEN
        DO I=1,NPARM
          CALL MNPOUT (I,CHNAM,VAL,ERROR,BND1,BND2,IVARBL)
          WRITE(22,1003)I,CHNAM,VAL,ERROR,BND1,BND2,IVARBL
 1003     FORMAT(7X,I5,1X,A10,4G20.10,I5)
        ENDDO
        CALL MNEMAT(EMAT,NPARM)
        DO I=1,NPARM
          WRITE(22,1002)(EMAT(I,J),J=1,NPARM)
 1002     FORMAT(8x,8G20.10)
        ENDDO
  
      ENDIF
      RETURN
      END        
      