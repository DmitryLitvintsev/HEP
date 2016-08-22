      PROGRAM MAIN
      COMMON /PAWC/ LPAW(300000)
      PARAMETER(MXBIN=51)
      COMMON /USER/ NDATA(MXBIN)
      EXTERNAL FCN
      PARAMETER(NCOLUMN=30)
      CHARACTER*5 CHTAGS(NCOLUMN)
      CHARACTER*80  CHTITL
      CHARACTER*132 HISNAM

      DATA CHTAGS /'VAR1', 'VAR2' ,'VAR3', 'VAR4', 'VAR5',
     +             'VAR6', 'VAR7' ,'VAR8', 'VAR9', 'VAR10',
     +             'VAR11','VAR12','VAR13','VAR14','VAR15',
     +             'VAR16','VAR17','VAR18','VAR19','VAR20',
     +             'VAR21','VAR22','VAR23','VAR24','VAR25',
     +             'VAR26','VAR27','VAR28','VAR29','VAR30'/
      DATA HISNAM /'misha.hist'/
     
      CALL HLIMIT(300000)
      
      CALL HROPEN(20,'paw',HISNAM,' ',1024,istat)
      IF (ISTAT.NE.0) WRITE(6,1001)HISNAM

C*
C*   Open data files
C*
      open(unit=5,file='fit.data',status='old')
      open(unit=6,file='fit.out',status='unknown',form='formatted')
      open(unit=21,file='fcn.dat',status='unknown',form='formatted')


      CALL HLDIR(' ',' ')
      CALL HRIN(0,999999,0)

      DO 100 I=1,51
      NDATA(I)=HI(11,I+13)
 100  CONTINUE

      CALL MINTIO(5,6,7)
      CALL MINUIT(FCN,0)
C*    CALL HROUT(0,ICYCLE,' ')
      CALL HREND('paw')
 1001 FORMAT(5X,'THE FILE NOT FOUND ',A132)
      STOP
      END

      SUBROUTINE FCN(NPAR,GIN,F,U,IFLAG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(*),GIN(*)
      PARAMETER (MXBIN=51)
      COMMON /USER/ NDATA(MXBIN)
C
C....FIT LMC2PI 'MISHA'S' CUTS BWIDTH=3MEV(HIST N 11 MISHA.HIST)
C
      DATA NBINS/51/
      IF (IFLAG.NE.1) GOTO 100
 100  SUM=0.D0
      DO 200 I=1,NBINS
         VAR=2.589D0+(DFLOAT(I)-0.5D0)*3.D-3
         FUN=BWG(VAR,U)
 200     SUM=SUM-PS(FUN,NDATA(I))
      F=SUM*2.D0
      IF (IFLAG.EQ.3) WRITE(21,1000  )F
1000  FORMAT( G20.10 )
      RETURN
      END
 
      FUNCTION FOLD(OBM,GAM,SIG,TRM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  FUNCTION FOLD - CONVOLUTION OF GAUSS AND BREIT-WIGNER
C    ESPESIALLY FOR PHI-MESON FITTING
C
C
C * INPUT: GAM,SIG,TRM,OBM;       OUTPUT: Y = FOLD(OBM,GAM,SIG,TRM);
C
C --- VARIABLES:
C          SIG : DETECTOR RESOLUTION
C          GAM : WIDTH OF THE RESONANCE
C          TRM : TRUE MASS
C          OBM : OBSERVED MASS
C          CWERF(Z) IS A PROGRAM OF THE CERN PROGRAM LIBRARY GENLIB
C
      COMPLEX*16 Z,WZ,WWERF
C------------------SUBSTITUTE
      X = (OBM - TRM) / (SIG * 1.414213D0)
      Y = GAM / (2.828427D0 * SIG)
      Z = CMPLX(X,Y)
C------------------COMPLEX ERF
      WZ = WWERF(Z)
C------------------GET NORM
      FOLD = DREAL(WZ)
      FOLD = FOLD / (2.50662827D0 * SIG )
      RETURN
      END
 
      FUNCTION BWG(X,U)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(*)
      SIG=     U(2)
      GAM=     U(3)
      TRM=U(1)
      A0=U(5)
      A1=U(6)
      A2=U(7)
      OBM=X
      BIN=3.D-3
      ATHR=2.589D0
      IF (X.LT.ATHR) THEN
         BWGP=0.D0
         RETURN
      ENDIF
 
      BWG=DABS(U(4)*BIN*FOLD(OBM,DABS(GAM),DABS(SIG)+1.E-20
     +                                           ,TRM))
     + +DABS(A0*(1.D0+A1*(X-TRM)
     + +A2*(X-TRM)**2))+1.D-31
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
