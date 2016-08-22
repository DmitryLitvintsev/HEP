HP FORTRAN 77                               Mon Jul 22 00:54:32 1996 
Ver: 09.16                                  ul.f                     Page     1

    1	      PROGRAM MAIN
    2	      COMMON /PAWC/ LPAW(300000)
    3	      PARAMETER(MXBIN=51)
    4	      COMMON /USER/ NDATA(MXBIN)
    5	      EXTERNAL FCN
    6	      PARAMETER(NCOLUMN=30)
    7	      CHARACTER*5 CHTAGS(NCOLUMN)
    8	      CHARACTER*80  CHTITL
    9	      CHARACTER*132 HISNAM
   10	
   11	      DATA CHTAGS /'VAR1', 'VAR2' ,'VAR3', 'VAR4', 'VAR5',
   12	     +             'VAR6', 'VAR7' ,'VAR8', 'VAR9', 'VAR10',
   13	     +             'VAR11','VAR12','VAR13','VAR14','VAR15',
   14	     +             'VAR16','VAR17','VAR18','VAR19','VAR20',
   15	     +             'VAR21','VAR22','VAR23','VAR24','VAR25',
   16	     +             'VAR26','VAR27','VAR28','VAR29','VAR30'/
   17	      DATA HISNAM /'misha.hist'/
   18	     
   19	      CALL HLIMIT(300000)
   20	      
   21	      CALL HROPEN(20,'paw',HISNAM,' ',1024,istat)
   22	      IF (ISTAT.NE.0) WRITE(6,1001)HISNAM
   23	
   24	C*
   25	C*   Open data files
   26	C*
   27	      open(unit=5,file='fit.data',status='old')
   28	      open(unit=6,file='fit.out',status='unknown',form='formatted')
   29	      open(unit=21,file='fcn.dat',status='unknown',form='formatted')
   30	
   31	
   32	      CALL HLDIR(' ',' ')
   33	      CALL HRIN(0,999999,0)
   34	
   35	      DO 100 I=1,51
   36	      NDATA(I)=HI(11,I+13)
   37	 100  CONTINUE
   38	
   39	      CALL MINTIO(5,6,7)
   40	      CALL MINUIT(FCN,0)
   41	C*    CALL HROUT(0,ICYCLE,' ')
   42	      CALL HREND('paw')
   43	 1001 FORMAT(5X,'THE FILE NOT FOUND ',A132)
   44	      STOP
   45	      END
**** WARNING # 1: on line 45 local variable chtitl never used
HP FORTRAN 77                               Mon Jul 22 00:54:32 1996 
Ver: 09.16                                  ul.f                     Page     2



Compilation statistics for procedure: main

Number of errors:       0		Number of Warnings:     1
Procedure number:     1

Accumulated number of source lines read:     45

   46	
   47	      SUBROUTINE FCN(NPAR,GIN,F,U,IFLAG)
   48	      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
   49	      DIMENSION U(*),GIN(*)
   50	      PARAMETER (MXBIN=51)
   51	      COMMON /USER/ NDATA(MXBIN)
   52	C
   53	C....FIT LMC2PI 'MISHA'S' CUTS BWIDTH=3MEV(HIST N 11 MISHA.HIST)
   54	C
   55	      DATA NBINS/51/
   56	      IF (IFLAG.NE.1) GOTO 100
   57	 100  SUM=0.D0
   58	      DO 200 I=1,NBINS
   59	         VAR=2.589D0+(DFLOAT(I)-0.5D0)*3.D-3
   60	         FUN=BWG(VAR,U)
   61	 200     SUM=SUM-PS(FUN,NDATA(I))
   62	      F=SUM*2.D0
   63	      IF (IFLAG.EQ.3) WRITE(21,1000  )F
   64	1000  FORMAT( G20.10 )
   65	      RETURN
   66	      END


Compilation statistics for procedure: fcn_

Number of errors:       0		Number of Warnings:     0
Procedure number:     2

Accumulated number of source lines read:     66

   67	 
   68	      REAL FUNCTION FOLD(OBM,GAM,SIG,TRM)
   69	C
   70	C  FUNCTION FOLD - CONVOLUTION OF GAUSS AND BREIT-WIGNER
   71	C    ESPESIALLY FOR PHI-MESON FITTING
   72	C
   73	C
   74	C * INPUT: GAM,SIG,TRM,OBM;       OUTPUT: Y = FOLD(OBM,GAM,SIG,TRM);
   75	C
   76	C --- VARIABLES:
   77	C          SIG : DETECTOR RESOLUTION
   78	C          GAM : WIDTH OF THE RESONANCE
   79	C          TRM : TRUE MASS
   80	C          OBM : OBSERVED MASS
HP FORTRAN 77                               Mon Jul 22 00:54:32 1996 
Ver: 09.16                                  ul.f                     Page     3

   81	C          CWERF(Z) IS A PROGRAM OF THE CERN PROGRAM LIBRARY GENLIB
   82	C
   83	      COMPLEX Z,WZ
   84	C------------------SUBSTITUTE
   85	      X = (OBM - TRM) / (SIG * 1.414213)
   86	      Y = GAM / (2.828427 * SIG)
   87	      Z = CMPLX(X,Y)
   88	C------------------COMPLEX ERF
   89	      WZ = CWERF(Z)
   90	C------------------GET NORM
   91	      FOLD = REAL(WZ)
   92	      FOLD = FOLD / (2.50662827 * SIG )
   93	      RETURN
   94	      END


Compilation statistics for procedure: fold_

Number of errors:       0		Number of Warnings:     0
Procedure number:     3

Accumulated number of source lines read:     94

   95	 
   96	      FUNCTION BWG(X,U)
   97	      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
   98	      DIMENSION U(*)
   99	      SIG=     U(2)
  100	      GAM=     U(3)
  101	      TRM=U(1)
  102	      A0=U(5)
  103	      A1=U(6)
  104	      A2=U(7)
  105	      OBM=X
  106	      BIN=3.D-3
  107	      ATHR=2.589D0
  108	      IF (X.LT.ATHR) THEN
  109	         BWGP=0.D0
  110	         RETURN
  111	      ENDIF
  112	 
  113	      BWG=DABS(U(4)*BIN*FOLD(OBM,DABS(GAM),DABS(SIG)+1.E-20
  114	     +                                           ,TRM))
  115	     + +DABS(A0*(1.D0+A1*(X-TRM)
  116	     + +A2*(X-TRM)**2))+1.D-31
  117	      END
HP FORTRAN 77                               Mon Jul 22 00:54:32 1996 
Ver: 09.16                                  ul.f                     Page     4



Compilation statistics for procedure: bwg_

Number of errors:       0		Number of Warnings:     0
Procedure number:     4

Accumulated number of source lines read:    117

  118	 
  119	      FUNCTION PS(X,N)
  120	      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  121	      IF (X.GT.0.D0 .AND. N.GT.1) THEN
  122	         S1 = DFLOAT(N) * DLOG(X)
  123	         S2 = 0.D0
  124	         DO 10 I = 2,N
  125	10       S2 = S2 + DLOG(DFLOAT(I))
  126	         PS = S1 - S2 - X
  127	      ENDIF
  128	      IF (X.EQ.0.D0             ) PS= 1.D0
  129	      IF (X.GT.0.D0 .AND. N.EQ.0) PS=-1.D0 * X
  130	      IF (X.GT.0.D0 .AND. N.EQ.1) PS=DLOG(X) - X
  131	      IF (X.LT.0                ) THEN
  132	          PRINT *,'******* FUNCTION LESS THAN 0 ********'
  133	          PRINT *,'******* PS IS SET TO -100    ********'
  134	          PS = -100.D0
  135	      ENDIF
  136	      RETURN
  137	      END


Compilation statistics for procedure: ps_

Number of errors:       0		Number of Warnings:     0
Procedure number:     5

Accumulated number of source lines read:    137



Final Compilation Statistics for file: ul.f

Accumulated errors:     0		Number of Warnings:     1
Total number of procedures:     5

Accumulated number of source lines read:    137


Timing Statistics for: 

f77pass1:	real:   0.35s   user:   0.00s    sys:   0.00s

