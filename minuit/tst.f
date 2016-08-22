      program main
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      external WWERF
      X = 10
      Y = 20 
      Z = CMPLX(X,Y)
      WZ = WWERF(Z)
      FOLD = DREAL(WZ)
      print*,'fold',fold 
      stop
      end
