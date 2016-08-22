#include "TH1.h"
#include "TF1.h"

void deconv(TH1* h,double hLow=-1, double hHigh=-1, double mass=1.321) {
	
	double hig = ( hHigh > 0 ) ? hHigh : h->GetXaxis()->GetXmax();
	double low = ( hLow > 0  ) ? hLow  : h->GetXaxis()->GetXmin();
	double bw  = h->GetBinWidth(0);
	
	TF1* fitFun = new TF1("deconv",deconv,low,hig,5);

  fitFun->SetLineWidth(2);
  fitFun->SetLineColor(4);
  
  	
  fitFun->SetParameter(0,0.5459);
  fitFun->SetParameter(1,-1.2427);
  fitFun->SetParameter(2,-2.9084);
  fitFun->SetParLimits(0,1,1);	
  fitFun->SetParLimits(1,1,1);	
  fitFun->SetParLimits(2,1,1);	
 		
  fitFun->SetParameter(3,1);
  fitFun->SetParameter(4,1);

  fitFun->SetRange(low,hig);

  h->Fit("deconv","QMR");
  h->Fit("deconv","QMR");
  h->Fit("deconv","EMQR");
  h->Fit("deconv","EMQR");
}
