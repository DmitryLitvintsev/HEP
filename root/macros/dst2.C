#include "TH1.h"
#include "TF1.h"

void dst2(TH1* h,int npar=3,double hLow=-1, double hHigh=-1, double mass=2.46,double sigma=0.020, char opt=' ') {

  double hig = ( hHigh > 0 ) ? hHigh : h->GetXaxis()->GetXmax();
  double low = ( hLow > 0  ) ? hLow  : h->GetXaxis()->GetXmin();
  double bw  = h->GetBinWidth(0);

  TF1* fitFun = new TF1("dst2",dst2,low,hig,4+npar);

  fitFun->SetLineWidth(2);
  fitFun->SetLineColor(4);
  
  
  fitFun->SetParameter(0,0);
  fitFun->SetParameter(1,bw);
  fitFun->SetParameter(2,sigma);
  fitFun->SetParameter(3,mass);

  for (int i=4; i<4+npar; i++) {
	  fitFun->SetParameter(i,0);
  }

  fitFun->SetRange(low,hig);

  fitFun->SetParLimits(0,1,1);
  fitFun->SetParLimits(1,1,1);
  fitFun->SetParLimits(2,1,1);
  fitFun->SetParLimits(3,1,1);
  
  fitFun->SetParName(0,"Entries");
  fitFun->SetParName(1,"binwidth");
  fitFun->SetParName(2,"sigma ");
  fitFun->SetParName(3,"Mass ");
  
  if (opt=='L') { 
	  h->Fit("dst2","QLMR");
	  fitFun->SetParLimits(0,0,1000000);
	  fitFun->SetParLimits(3,mass-10.*sigma,mass+10.*sigma);
	  h->Fit("dst2","QLMR");
	  fitFun->SetParLimits(2,sigma/10.,sigma*5.);
	  h->Fit("dst2","QLMR");
	  h->Fit("dst2","ELMQR");
	  
  }
  else { 
	  h->Fit("dst2","QMR");
	  fitFun->SetParLimits(0,0,1000000);
	  fitFun->SetParLimits(3,mass-10.*sigma,mass+10.*sigma);
	  h->Fit("dst2","QMR");
	  fitFun->SetParLimits(2,sigma/10.,sigma*5.);
	  h->Fit("dst2","QMR");
	  h->Fit("dst2","EMVR");
  }
  delete fitFun;
}
