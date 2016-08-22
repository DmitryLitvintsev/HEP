#include "TH1.h"
#include "TF1.h"

void effic(TH1* h,double hLow=-1, double hHigh=-1, double mass=1.321) {

	double hig = ( hHigh > 0 ) ? hHigh : h->GetXaxis()->GetXmax();
	double low = ( hLow > 0  ) ? hLow  : h->GetXaxis()->GetXmin();
	double bw  = h->GetBinWidth(0);

	TF1* fitFun = new TF1("eff",eff,low,hig,2);

  fitFun->SetLineWidth(2);
  fitFun->SetLineColor(4);
  
  
  fitFun->SetParameter(0,1);
  fitFun->SetParameter(1,1);
//  fitFun->SetParameter(2,1);

  fitFun->SetRange(low,hig);

  h->Fit("eff","QMR");
  h->Fit("eff","QMR");
  h->Fit("eff","EMQR");
  h->Fit("eff","EMVR");
}
