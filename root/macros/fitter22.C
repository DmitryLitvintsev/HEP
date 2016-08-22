#include "TH1.h"
#include "TF1.h"

void fitter22(TH1* h,
	      int npar=3,
	      double hLow=-1, 
	      double hHigh=-1, 
	      char ch=' ') {

	double hLow=-1;
	double hHigh=-1;

	int npar=3;  
  double hig = ( hHigh > 0 ) ? hHigh : h->GetXaxis()->GetXmax();
  double low = ( hLow > 0  ) ? hLow  : h->GetXaxis()->GetXmin();

  double bw  = h->GetBinWidth(0);

  TF1* fitFun = new TF1("gp2",gp22,low,hig,7+npar);

{  for (int i=7; i<7+npar; i++) {
	  fitFun->SetParameter(i,0);
  }
}


  fitFun->SetLineWidth(2);
  fitFun->SetLineColor(4);
  
  double sigma=0.02;
  
  fitFun->SetParameter(0,0);
  fitFun->SetParameter(1,bw);
  fitFun->SetParameter(2,sigma);
  fitFun->SetParameter(3,5.792);
  fitFun->SetParameter(4,0);
  fitFun->SetParameter(5,0.595);
  fitFun->SetParameter(6,4.25);

  fitFun->SetRange(low,hig);

  fitFun->SetParLimits(0,1,1);
  fitFun->SetParLimits(1,1,1);
  fitFun->SetParLimits(2,1,1);
  fitFun->SetParLimits(3,1,1);
  fitFun->SetParLimits(4,1,1);
  fitFun->SetParLimits(5,1,1);
  
  fitFun->SetParName(0,"n1");
  fitFun->SetParName(1,"binwidth");
  fitFun->SetParName(2,"sigma1 ");
  fitFun->SetParName(5,"Mass1");
  fitFun->SetParName(4,"n2");
  fitFun->SetParName(5,"sigma2");
  fitFun->SetParName(6,"Mass2");
  fitFun->SetParName(4,"n2");

  if ( ch!='L') { 
	  h->Fit("gp2","QMR");
	  h->Fit("gp2","QMR");
	  fitFun->ReleaseParameter(0);
	  fitFun->ReleaseParameter(3);
 	  fitFun->SetParLimits(5,mass-10.*sigma,mass+10.*sigma);
	  h->Fit("gp2","QMR");
	  h->Fit("gp2","QMR");
 	  fitFun->SetParLimits(2,sigma/10.,sigma*2.);
 	  fitFun->SetParLimits(4,sigma/10.,sigma*2.);
	  h->Fit("gp2","MQR");
	  h->Fit("gp2","EMR");
// 	  h->Fit("gp2","EMQR");
// 	  h->Fit("gp2","EMQR");
// 	  h->Fit("gp2","EMQR");
// 	  h->Fit("gp2","EMQR");  
  }
  else { 
	  h->Fit("gp2","LQMR");
	  h->Fit("gp2","LQMR");
	  fitFun->SetParLimits(0,0,1000000);
	  fitFun->SetParLimits(3,0,1000000);
 	  fitFun->SetParLimits(5,mass-10.*sigma,mass+10.*sigma);
	  h->Fit("gp2","LQMR");
	  h->Fit("gp2","LQMR");
 	  fitFun->SetParLimits(2,sigma/10.,sigma*5.);
 	  fitFun->SetParLimits(4,sigma/10.,sigma*5.);
	  h->Fit("gp2","LMQR");
	  h->Fit("gp2","LEMR");
// 	  h->Fit("gp2","LEMQR");
// 	  h->Fit("gp2","LEMQR");
// 	  h->Fit("gp2","LEMQR");
// 	  h->Fit("gp2","LEMQR");
  }
  delete fitFun;
}
