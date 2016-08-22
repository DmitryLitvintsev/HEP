      PROGRAM MAIN
      IMPLICIT NONE
      INTEGER LPAW
      COMMON /PAWC/ LPAW(300000)
      DOUBLE PRECISION U, LH
      DIMENSION U(*)
      

      CALL HLIMIT(300000)
      
      CALL HROPEN(20,'blah','lh.hist','N',1024,istat)
      CALL HBOOK1(1,'xicst',700.0,1.4,2.8,0.0)
c
      
