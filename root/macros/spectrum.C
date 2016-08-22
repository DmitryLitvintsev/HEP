#include "TH1.h"
#include "TF1.h"

void spectrum(TH1* h,double hLow=-1, double hHigh=-1, char opt=' ') {
	
	double hig = ( hHigh > 0 ) ? hHigh : h->GetXaxis()->GetXmax();
	double low = ( hLow > 0  ) ? hLow  : h->GetXaxis()->GetXmin();
	double bw  = h->GetBinWidth(0);

	
	TF1* fitFun = new TF1("eff",GauExp,low,hig,4);

	fitFun->SetLineWidth(2);
	fitFun->SetLineColor(4);
  
	
	fitFun->SetParameter(0,1);
	fitFun->SetParameter(1,0.86);
	fitFun->SetParameter(2,3.8);
	fitFun->SetParameter(3,1.);
	fitFun->SetRange(low,hig);
	
	if (opt=='L') { 
		h->Fit("eff","LQMR");
		h->Fit("eff","LQMR");
		h->Fit("eff","LEMQR");
		h->Fit("eff","LEVMR");
	}
	else { 
		h->Fit("eff","QMR");
		h->Fit("eff","QMR");
		h->Fit("eff","EMQR");
		h->Fit("eff","EVMR");
	}
}
