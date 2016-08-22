#include "TH1.h"
#include "TF1.h"
void ddsbkg(TH1* h,double hLow=-1, double hHigh=-1,char tag=' ') {

	double hig = ( hHigh > 0 ) ? hHigh : h->GetXaxis()->GetXmax();
	double low = ( hLow > 0  ) ? hLow  : h->GetXaxis()->GetXmin();
	double bw  = h->GetBinWidth(0);

	TF1* fitFun = new TF1("ddstbkg",ddst,low,hig,4);

	fitFun->SetLineWidth(2);
	fitFun->SetLineColor(4);
  
	if (tag == 'L') { 
		h->Fit("ddst","LMR");
		h->Fit("ddst","ELMR");
	}
	else { 
		h->Fit("ddst","MR");
		h->Fit("ddst","EMR");
	}
	
	delete fitFun;
}
