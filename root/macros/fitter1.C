#include "TH1.h"
#include "TF1.h"

void fitter1(TH1* h,
	     int npar=3,
	     double hLow=-1, 
	     double hHigh=-1, 
	     double mass=2.46,
	     double sigma=0.020, 
	     char opt=' ',
	     int number=0  ) {

//  double hig = ( hHigh > 0 ) ? hHigh : h->GetXaxis()->GetXmax();
//  double low = ( hLow > 0  ) ? hLow  : h->GetXaxis()->GetXmin();
  double hig =  hHigh;
  double low =  hLow;
  double bw  = h->GetBinWidth(0);

  TF1* fitFun = new TF1("gp",gp1,low,hig,4+npar);

  fitFun->SetLineWidth(2);
  fitFun->SetLineColor(4);
  
  
  fitFun->SetParameter(0,number);
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
 	  h->Fit("gp","LMR");
 	  fitFun->ReleaseParameter(0);
 	  fitFun->ReleaseParameter(3);
 	  h->Fit("gp","LMR");
 	  fitFun->ReleaseParameter(2);
 	  h->Fit("gp","LEMR");
   }
   else { 
   	  h->Fit("gp","MR");
   	  fitFun->ReleaseParameter(0);
   	  fitFun->ReleaseParameter(3);
   	  h->Fit("gp","MR");
   	  fitFun->ReleaseParameter(2);
   	  h->Fit("gp","MR");
   	  h->Fit("gp","EMR");
   }

  delete fitFun;
}
