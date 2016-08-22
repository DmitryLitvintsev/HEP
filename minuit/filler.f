      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER LPAW,ICYCLE,ISTAT
      COMMON /PAWC/ LPAW(300000)
      REAL Y,X,EX,LH
      DIMENSION LH(2000)
      INTEGER COUNT,RUN,EVENT
      CALL HLIMIT(300000)
      CALL HROPEN(20,'blah','tst1.hist','N',1024,ISTAT)

      CALL HBOOK1(1,'data',280.0,1.4,2.8,0.0)

C      CALL HBOOK1(1,'LH',2000,-1000.0,1000.0,0.0)
C      CALL HBOOK1(2,'profile',2000,-1000.0,1000.0,0.0)

C      Y = 0.0

C      OPEN(UNIT=80,FILE='fcn.dat')
C      DO I=-1000,1000
C        READ(80,*,END=500)Y
C        IF (FLOAT(I).GE.0) THEN
C        ENDIF
C      ENDDO
C 500  CONTINUE

      OPEN(UNIT=82,FILE='data_points_right_sign.dat')
 52   READ(82,*,END=502)COUNT,RUN,EVENT,X,EX
      CALL HF1(1,X,1.0)
      GOTO 52
 502  CONTINUE



      CALL HROUT (0,ICYCLE,' ')
      CALL HREND('blah')
      STOP
      END
            
