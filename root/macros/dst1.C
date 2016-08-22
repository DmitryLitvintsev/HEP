#include "TH1.h"
#include "TF1.h"

void dst1(TH1* h,int npar=3,double hLow=-1, double hHigh=-1, double mass=1.321,double sigma=0.003,char ch=' ') {

  double hig = ( hHigh > 0 ) ? hHigh : h->GetXaxis()->GetXmax();
  double low = ( hLow > 0  ) ? hLow  : h->GetXaxis()->GetXmin();
  double bw  = h->GetBinWidth(0);

  TF1* fitFun = new TF1("dstfit",dstfit,low,hig,6+npar);

  for (int i=6; i<6+npar; i++) {
	  fitFun->SetParameter(i,0);
  }


  fitFun->SetLineWidth(2);
  fitFun->SetLineColor(4);
  
  
  fitFun->SetParameter(0,0);
  fitFun->SetParameter(1,bw);
  fitFun->SetParameter(2,sigma);
  fitFun->SetParameter(3,0);
  fitFun->SetParameter(4,sigma);
  fitFun->SetParameter(5,mass);

//    fitFun->SetParameter(3,1.9);
//    fitFun->SetParameter(2,0.0036);
//    fitFun->SetParameter(4,0.0016);


//    fitFun->SetParameter(3,1.2);
//    fitFun->SetParameter(2,0.0061);
//    fitFun->SetParameter(4,0.0023);


  fitFun->SetRange(low,hig);

  fitFun->SetParLimits(0,1,1);
  fitFun->SetParLimits(1,1,1);
  fitFun->SetParLimits(2,1,1);
  fitFun->SetParLimits(3,1,1);
  fitFun->SetParLimits(4,1,1);
  fitFun->SetParLimits(5,1,1);
  
  fitFun->SetParName(0,"n1");
  fitFun->SetParName(3,"ratio");
  fitFun->SetParName(1,"binwidth");
  fitFun->SetParName(2,"sigma1 ");
  fitFun->SetParName(5,"Mass1");
  fitFun->SetParName(4,"sigma2");

  if ( ch!='L') { 
	  h->Fit("dstfit","MR");
	  fitFun->SetParLimits(0,0,1000000);
	  fitFun->SetParLimits(3,0,1000000);
	  fitFun->SetParLimits(5,mass-10.*sigma,mass+10.*sigma);
	  h->Fit("dstfit","MR");
	  fitFun->SetParLimits(2,sigma/100.,sigma*10.);
	  fitFun->SetParLimits(4,sigma/100.,sigma*10.);
	  h->Fit("dstfit","MR");
	  h->Fit("dstfit","EMQR");
// 	  h->Fit("dstfit","EMQR");
// 	  h->Fit("dstfit","EMQR");
// 	  h->Fit("dstfit","EMQR");
// 	  h->Fit("dstfit","EMQR");
  }
  else { 
	  h->Fit("dstfit","LQMR");
	  h->Fit("dstfit","LQMR");
	  fitFun->SetParLimits(0,0,1000000);
	  fitFun->SetParLimits(3,0,1000000);
 	  fitFun->SetParLimits(5,mass-10.*sigma,mass+10.*sigma);
	  h->Fit("dstfit","LQMR");
	  h->Fit("dstfit","LQMR");
 	  fitFun->SetParLimits(2,sigma/10.,sigma*2.);
 	  fitFun->SetParLimits(4,sigma/10.,sigma*2.);
	  h->Fit("dstfit","LMQR");
	  h->Fit("dstfit","LEMQR");
// 	  h->Fit("dstfit","LEMQR");
// 	  h->Fit("dstfit","LEMQR");
// 	  h->Fit("dstfit","LEMQR");
// 	  h->Fit("dstfit","LEMQR");
  }
  delete fitFun;
}
