#include "TH1.h"
#include "TF1.h"

void penta2(TH1F* h,
	   int npar=3,
	   double hLow=-1, 
	   double hHigh=-1, 
	   double mass=1.53, 
	   double sigma=0.003, 
	   int number=0) {

  double hig = ( hHigh > 0 ) ? hHigh : h->GetXaxis()->GetXmax();
  double low = ( hLow > 0  ) ? hLow  : h->GetXaxis()->GetXmin();
  double bw  = h->GetBinWidth(0);

   TF1* fitFun = new TF1("pentafit1",pentafit1,low,hig,8+npar);

  for (int i=8; i<8+npar; i++) {
	  fitFun->SetParameter(i,0);
  }


  fitFun->SetLineWidth(2);
  fitFun->SetLineColor(4);
  
  fitFun->SetParameter(0,2168);
  fitFun->SetParameter(1,bw);
  fitFun->SetParameter(2,0.0057);
  fitFun->SetParameter(3,0.0091);
  fitFun->SetParameter(4,1.5320);
  fitFun->SetParameter(5,number);
  fitFun->SetParameter(6,0.008);
  fitFun->SetParameter(7,1.862);
  fitFun->SetParameter(8,1265.7126);
  fitFun->SetParameter(9,-1.5855);
  fitFun->SetParameter(10,2.0099);
  fitFun->SetParameter(11,-0.5594);

  fitFun->SetParName(0,"n1");
  fitFun->SetParName(1,"binwidth");
  fitFun->SetParName(2,"sigma1 ");
  fitFun->SetParName(3,"Gamma ");
  fitFun->SetParName(4,"M1 ");
  fitFun->SetParName(5,"n2");
  fitFun->SetParName(6,"sigma2 ");
  fitFun->SetParName(7,"M2 ");

  fitFun->SetRange(low,hig);

//  fitFun->SetParLimits(0,1,1);
  fitFun->SetParLimits(1,1,1);
//  fitFun->SetParLimits(2,1,1);
  fitFun->SetParLimits(3,1,1);
//  fitFun->SetParLimits(4,1,1);
//  fitFun->SetParLimits(5,1,1);
  fitFun->SetParLimits(6,1,1);
  fitFun->SetParLimits(7,1,1);
  fitFun->SetParLimits(8,1,1);
  fitFun->SetParLimits(9,1,1);
  fitFun->SetParLimits(10,1,1);
  h->Fit("pentafit1","QMR");
  h->Fit("pentafit1","EMQR");
  delete fitFun;
  
}

