#include "TH1.h"
#include "TF1.h"

void dstexp(TH1* h,double hLow=-1, double hHigh=-1, double mass=1.321,double sigma=0.003,char ch=' ') {

  double hig = ( hHigh > 0 ) ? hHigh : h->GetXaxis()->GetXmax();
  double low = ( hLow > 0  ) ? hLow  : h->GetXaxis()->GetXmin();
  double bw  = h->GetBinWidth(0);

  TF1* fitFun = new TF1("dstfitexp",dstfitexp,low,hig,10);


  fitFun->SetLineWidth(2);
  fitFun->SetLineColor(4);
  
  
  fitFun->SetParameter(0,0);
  fitFun->SetParameter(1,bw);
  fitFun->SetParameter(2,sigma);
  fitFun->SetParameter(3,0);
  fitFun->SetParameter(4,sigma);
  fitFun->SetParameter(5,mass);
  fitFun->SetParameter(6,mass);
  fitFun->SetParameter(7,1.);
  fitFun->SetParameter(8,1.);
  fitFun->SetParameter(9,2.);

  

  
  fitFun->SetRange(low,hig);

  fitFun->SetParLimits(0,1,1);
  fitFun->SetParLimits(1,1,1);
  fitFun->SetParLimits(2,1,1);
  fitFun->SetParLimits(3,1,1);
  fitFun->SetParLimits(4,1,1);
  fitFun->SetParLimits(5,1,1);
  fitFun->SetParLimits(6,1,1);
  
  fitFun->SetParName(0,"n1");
  fitFun->SetParName(3,"ratio");
  fitFun->SetParName(1,"binwidth");
  fitFun->SetParName(2,"sigma1 ");
  fitFun->SetParName(5,"Mass1");
  fitFun->SetParName(6,"Mass2");
  fitFun->SetParName(4,"sigma2");

  if ( ch!='L') { 
	  h->Fit("dstfitexp","MR");
	  fitFun->SetParLimits(0,0,1000000);
	  fitFun->SetParLimits(3,0,1000000);
	  fitFun->SetParLimits(5,mass-10.*sigma,mass+10.*sigma);
	  fitFun->SetParLimits(6,mass-10.*sigma,mass+10.*sigma);
	  h->Fit("dstfitexp","MR");
	  fitFun->SetParLimits(2,sigma/100.,sigma*10.);
	  fitFun->SetParLimits(4,sigma/100.,sigma*10.);
	  h->Fit("dstfitexp","MR");
	  h->Fit("dstfitexp","EMQR");
// 	  h->Fit("dstfitexp","EMQR");
// 	  h->Fit("dstfitexp","EMQR");
// 	  h->Fit("dstfitexp","EMQR");
// 	  h->Fit("dstfitexp","EMQR");
  }
  else { 
	  h->Fit("dstfitexp","LQMR");
	  h->Fit("dstfitexp","LQMR");
	  fitFun->SetParLimits(0,0,1000000);
	  fitFun->SetParLimits(3,0,1000000);
 	  fitFun->SetParLimits(5,mass-10.*sigma,mass+10.*sigma);
 	  fitFun->SetParLimits(6,mass-10.*sigma,mass+10.*sigma);
	  h->Fit("dstfitexp","LQMR");
	  h->Fit("dstfitexp","LQMR");
 	  fitFun->SetParLimits(2,sigma/10.,sigma*2.);
 	  fitFun->SetParLimits(4,sigma/10.,sigma*2.);
	  h->Fit("dstfitexp","LMQR");
	  h->Fit("dstfitexp","LEMQR");
// 	  h->Fit("dstfitexp","LEMQR");
// 	  h->Fit("dstfitexp","LEMQR");
// 	  h->Fit("dstfitexp","LEMQR");
// 	  h->Fit("dstfitexp","LEMQR");
  }
  delete fitFun;
}
