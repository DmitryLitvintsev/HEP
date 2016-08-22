#include "TH1.h"
#include "TF1.h"

void sigmacfitter(TH1* h,int npar=3,double hLow=-1, double hHigh=-1, double mass=1.321,double sigma=0.003,char ch=' ') {

  double hig = ( hHigh > 0 ) ? hHigh : h->GetXaxis()->GetXmax();
  double low = ( hLow > 0  ) ? hLow  : h->GetXaxis()->GetXmin();
  double bw  = h->GetBinWidth(0);

  TF1* fitFun = h->GetFunction("sigmacfit");

  if (!fitFun) { 
  
	  fitFun = new TF1("sigmacfit",sigmacfit,low,hig,9);
 

	  fitFun->SetLineWidth(2);
	  fitFun->SetLineColor(4);
	  
	  
	  fitFun->SetParameter(0,0);
	  fitFun->SetParameter(1,bw);
	  fitFun->SetParameter(2,1.25e-3);
	  fitFun->SetParameter(3,0.00);
	  fitFun->SetParameter(4,0.169);
	  fitFun->SetParameter(5,0.5);
	  fitFun->SetParameter(6,0.34);

	  fitFun->SetRange(low,hig);
	  
	  fitFun->SetParLimits(0,1,1);
	  fitFun->SetParLimits(1,1,1);
	  fitFun->SetParLimits(2,1,1);
	  fitFun->SetParLimits(3,1,1);
	  fitFun->SetParLimits(4,1,1);
	  fitFun->SetParLimits(5,1,1);
	  fitFun->SetParLimits(6,1,1);
//	  fitFun->SetParameter(2,1.25594e-03);
	  
	  fitFun->SetParName(0,"N1");
	  fitFun->SetParName(1,"binwidth");
	  fitFun->SetParName(2,"sigma");
	  fitFun->SetParName(3,"Gamma");
	  fitFun->SetParName(4,"M1");
	  fitFun->SetParName(5,"pwr");
	  
	  h->Fit("sigmacfit","QMR");
	  fitFun->ReleaseParameter(5);
	  h->Fit("sigmacfit","QMR");
	  fitFun->ReleaseParameter(6);
	  h->Fit("sigmacfit","QMR");
	  fitFun->ReleaseParameter(2); 
	  h->Fit("sigmacfit","EVMR");
	  fitFun->ReleaseParameter(0); 
//	  fitFun->ReleaseParameter(3); 
	  fitFun->ReleaseParameter(4); 
//	  fitFun->SetParameter(3,0.0022);
          h->Fit("sigmacfit","EVMR");
  }
  else { 
	  fitFun->ReleaseParameter(0); 
	  fitFun->ReleaseParameter(3); 
	  fitFun->SetParameter(3,0);	
	  fitFun->SetParLimits(3,1,1);
//	  fitFun->SetParLimits(2,1,1);
	  h->Fit("sigmacfit","EMR");
	  fitFun->ReleaseParameter(4); 
	  h->Fit("sigmacfit","EVMR");
  }
  
}
