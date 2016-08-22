#include "TH1.h"
#include "TF1.h"


//
// this macro was used to produce LH scan for (Xi-,pi-) G + SQRT(x)*Pol
//

void penta(TH1* h,
	   int npar=3,
	   double hLow=-1, 
	   double hHigh=-1, 
	   double mass=1.53, 
	   double sigma=0.003, 
	   int number=0) {

  double hig = ( hHigh > 0 ) ? hHigh : h->GetXaxis()->GetXmax();
  double low = ( hLow > 0  ) ? hLow  : h->GetXaxis()->GetXmin();
  double bw  = h->GetBinWidth(0);

   TF1* fitFun = new TF1("pentafit",pentafit,low,hig,4+npar);

  for (int i=4; i<4+npar; i++) {
	  fitFun->SetParameter(i,0);
  }


  fitFun->SetLineWidth(2);
  fitFun->SetLineColor(4);
  
  fitFun->SetParameter(0,number);
  fitFun->SetParameter(1,bw);
  fitFun->SetParameter(2,sigma);
  fitFun->SetParameter(3,mass);

  fitFun->SetParName(0,"Entries");
  fitFun->SetParName(1,"binwidth");
  fitFun->SetParName(2,"sigma ");
  fitFun->SetParName(3,"Mass ");

  fitFun->SetRange(low,hig);

  

  fitFun->SetParLimits(0,1,1);
  fitFun->SetParLimits(1,1,1);
  fitFun->SetParLimits(2,1,1);
  fitFun->SetParLimits(3,1,1);

  h->Fit("pentafit","QMR");
  fitFun->SetParLimits(0,-10000,1000000);
//  fitFun->SetParLimits(3,1.862-0.002,1.862+0.002);
  fitFun->SetParLimits(3,mass-10.*sigma,mass+10.*sigma);
  h->Fit("pentafit","QMR");
  fitFun->SetParLimits(2,sigma/10.,sigma*5.);
//  fitFun->SetParLimits(2,sigma-0.0012,sigma+0.0012);
  h->Fit("pentafit","QMR");
  h->Fit("pentafit","EQMR");
  delete fitFun;

 }
