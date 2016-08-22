#include "TH1.h"
#include "TF1.h"
void xicst(TH1* h,int npar=3,double hLow=-1, double hHigh=-1, double mass=1.53,double sigma=0.003) {

  double hig = ( hHigh > 0 ) ? hHigh : h->GetXaxis()->GetXmax();
  double low = ( hLow > 0  ) ? hLow  : h->GetXaxis()->GetXmin();
  double bw  = h->GetBinWidth(0);

  TF1* fitFun = new TF1("gbwp",xistfit,low,hig,9);

  fitFun->SetLineWidth(2);
  fitFun->SetLineColor(4);
  
  
  fitFun->SetParameter(0,0);
  fitFun->SetParameter(1,bw);
  fitFun->SetParameter(2,0.0036);
  fitFun->SetParameter(3,0.0091);
  fitFun->SetParameter(4,mass);
  fitFun->SetParameter(5,0);
  fitFun->SetParameter(6,0);

  fitFun->SetParameter(7,0);
  fitFun->SetParameter(8,0);
//   fitFun->FixParameter(7,0);
//   fitFun->FixParameter(8,0);

  fitFun->SetRange(low,hig);

  fitFun->SetParLimits(0,1,1);
  fitFun->SetParLimits(1,1,1);
  fitFun->SetParLimits(2,1,1);
  fitFun->SetParLimits(3,1,1);
  fitFun->SetParLimits(4,1,1);
  
  fitFun->SetParName(0,"Entries");
  fitFun->SetParName(1,"binwidth");
  fitFun->SetParName(2,"sigma ");
  fitFun->SetParName(4,"Mass ");
  fitFun->SetParName(3,"Gamma ");

  h->Fit("gbwp","QR");
  h->Fit("gbwp","QR");
  fitFun->SetParLimits(0,0,1000000);
  fitFun->SetParLimits(4,1.52,1.55);
  h->Fit("gbwp","QR");
  h->Fit("gbwp","QR");
//  fitFun->SetParLimits(2,0.0001,0.010);
  h->Fit("gbwp","EMVR");
  h->Fit("gbwp","EMVR");
  delete fitFun;
}
