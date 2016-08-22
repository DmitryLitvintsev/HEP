#include "TH1.h"
#include "TF1.h"
void lambda(TH1* h,int npar=3,double hLow=-1, double hHigh=-1, double mass=1.53,double sigma=0.003, char opt=' ') {

  double hig = ( hHigh > 0 ) ? hHigh : h->GetXaxis()->GetXmax();
  double low = ( hLow > 0  ) ? hLow  : h->GetXaxis()->GetXmin();
  double bw  = h->GetBinWidth(0);

  TF1* fitFun = new TF1("lambdafit",lambdafit,low,hig,6+npar);

  for (int i=6; i<6+npar; i++) {
	  fitFun->SetParameter(i,0);
  }

//  FCN=7846 FROM MINOS     STATUS=SUCCESSFUL    510 CALLS        3341 TOTAL
//                      EDM=6.71109e-08    STRATEGY= 1  ERROR MATRIX UNCERTAINTY   2.0 per cent
//   EXT PARAMETER                                   STEP         FIRST   
//   NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
//    1  Entries      1.32513e+06   1.89453e+03  -8.85395e-08   1.55867e-01
//    2  binwidth     5.00000e-04     fixed    
//    3  sigma1       1.20972e-03   4.63686e-06   1.35061e-09   2.27265e-01
//    4  p3           8.74640e-01   8.57157e-03  -5.47007e-14   4.16129e+05
//    5  sigma2       3.64497e-03   1.57722e-05  -6.67236e-08  -8.16851e-01
//    6  Mass1        1.11561e+00   2.14554e-06  -3.80138e-08   9.02944e-02
//    7  p6           7.09129e+04   1.58426e+02   4.82957e-02   2.24535e-07
//    8  p7           7.32996e+00   4.74657e-02   6.19256e-05   4.95304e-04
//    9  p8          -3.14016e+02   3.76458e+00   3.76458e+00   3.12905e-06

  fitFun->SetLineWidth(2);
  fitFun->SetLineColor(4);
  
  
  fitFun->SetParameter(0,0);
  fitFun->SetParameter(1,bw);
  fitFun->SetParameter(2,1.20972e-03);
  fitFun->SetParameter(4,3.64497e-03);
  fitFun->SetParameter(5,mass);
  fitFun->SetParameter(3,1);

  fitFun->SetRange(low,hig);

  fitFun->SetParLimits(0,1,1);
  fitFun->SetParLimits(1,1,1);
  fitFun->SetParLimits(2,1,1);
  fitFun->SetParLimits(3,1,1);
  fitFun->SetParLimits(4,1,1);
  fitFun->SetParLimits(5,1,1);
  
  fitFun->SetParName(0,"Entries");
  fitFun->SetParName(1,"binwidth");
  fitFun->SetParName(2,"sigma1");
  fitFun->SetParName(4,"sigma2");
  fitFun->SetParName(5,"Mass1");

  if (opt!='L') { 
	  h->Fit("lambdafit","QR");
	  fitFun->SetParLimits(0,0,100000000);
	  fitFun->SetParLimits(3,-100,100);
	  fitFun->SetParLimits(5,mass-10.*sigma,mass+10.*sigma);
	  h->Fit("lambdafit","VMR");
	  fitFun->SetParLimits(2,1.e-5,10);
	  fitFun->SetParLimits(4,1.e-5,10);
	  h->Fit("lambdafit","VEMR");
  }
  else { 
	  h->Fit("lambdafit","LQR");
	  fitFun->SetParLimits(0,0,100000000);
	  fitFun->SetParLimits(3,-100000000,100000000);
	  fitFun->SetParLimits(5,mass-10.*sigma,mass+10.*sigma);
	  h->Fit("lambdafit","LMVR");
	  fitFun->SetParLimits(2,1.e-5,10);
	  fitFun->SetParLimits(4,1.e-5,10);
	  h->Fit("lambdafit","LMVR");
	  h->Fit("lambdafit","LVEMR");
  }
  delete fitFun;
}
