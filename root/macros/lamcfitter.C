#include "TH1.h"
#include "TF1.h"

void lamcfitter(TH1* h,int npar=3,double hLow=-1, double hHigh=-1, double mass=1.321,double sigma=0.003,char ch=' ') {

  double hig = ( hHigh > 0 ) ? hHigh : h->GetXaxis()->GetXmax();
  double low = ( hLow > 0  ) ? hLow  : h->GetXaxis()->GetXmin();
  double bw  = h->GetBinWidth(0);

  TF1* fitFun = h->GetFunction("lamcstfit");
  if (fitFun==NULL) { 
  
	  fitFun = new TF1("lamcstfit",lamcstfit,low,hig,10);
 

	  fitFun->SetLineWidth(2);
	  fitFun->SetLineColor(4);
	  
	  
	  fitFun->SetParameter(0,0);
	  fitFun->SetParameter(1,bw);
	  fitFun->SetParameter(2,sigma);
	  fitFun->SetParameter(3,0.00);
	  fitFun->SetParameter(4,3.07999999999999829e-01);
	  fitFun->SetParameter(5,0);
	  fitFun->SetParameter(6,0.34);
	  fitFun->SetParameter(7,2.);  
	  fitFun->SetParameter(8,0.);
          fitFun->SetParameter(9,0.);
          fitFun->FixParameter(9,0.);
	  

	  fitFun->SetRange(low,hig);
	  
	  fitFun->SetParLimits(0,1,1);
	  fitFun->SetParLimits(1,1,1);
	  fitFun->SetParLimits(2,1,1);
	  fitFun->SetParLimits(3,1,1);
	  fitFun->SetParLimits(4,1,1);
	  fitFun->SetParLimits(5,1,1);
	  fitFun->SetParLimits(6,1,1);
	  
	  fitFun->SetParName(0,"N1");
	  fitFun->SetParName(1,"binwidth");
	  fitFun->SetParName(2,"sigma");
	  fitFun->SetParName(3,"Gamma");
	  fitFun->SetParName(4,"M1");
	  fitFun->SetParName(5,"N2");
	  fitFun->SetParName(6,"M2");
	  fitFun->SetParName(7,"pwr");
	  
	  h->Fit("lamcstfit","LQMR");
	  fitFun->ReleaseParameter(5);
	  h->Fit("lamcstfit","LQMR");
	  fitFun->ReleaseParameter(6);
	  h->Fit("lamcstfit","LQMR");
	  fitFun->ReleaseParameter(2); 
	  h->Fit("lamcstfit","LEVMR");
	  fitFun->ReleaseParameter(0); 
	  fitFun->ReleaseParameter(3); 
	  fitFun->ReleaseParameter(4); 
	  fitFun->SetParLimits(4,0.33,0.35);
	  fitFun->SetParameter(3,0.003);
//  h->Fit("lamcstfit","LEMR");
  }
  else { 
	  fitFun->ReleaseParameter(0); 
//	  fitFun->ReleaseParameter(3); 
	  fitFun->SetParameter(3,0.003);
	  h->Fit("lamcstfit","LEMR");
	  fitFun->ReleaseParameter(4); 
	  h->Fit("lamcstfit","LEMR");
  }
	  
}
