#include "TH1.h"
#include "TF1.h"
void ddst1(TH1* h,int npar=3,double hLow=-1, double hHigh=-1, double mass=1.53,double sigma=0.003) {

  double hig = ( hHigh > 0 ) ? hHigh : h->GetXaxis()->GetXmax();
  double low = ( hLow > 0  ) ? hLow  : h->GetXaxis()->GetXmin();
  double bw  = h->GetBinWidth(0);

  TF1* fitFun = new TF1("ddst",ddst,low,hig,11);

  fitFun->SetLineWidth(2);
  fitFun->SetLineColor(4);
  
  
  fitFun->SetParameter(0,0);
  fitFun->SetParameter(1,bw);
  fitFun->SetParameter(2,0.0038);
  fitFun->SetParameter(3,0.020);
  fitFun->SetParameter(4,0.413);

  fitFun->SetParameter(9,0.0);
  fitFun->SetParameter(10,0.020);
  fitFun->SetParameter(11,0.455);
  fitFun->SetParameter(4,mass);

  fitFun->SetParameter(5,0);
  fitFun->SetParameter(6,0);
  fitFun->SetParameter(7,0);
  fitFun->SetParameter(8,0);

  fitFun->SetRange(low,hig);

  fitFun->SetParLimits(0,1,1);
  fitFun->SetParLimits(1,1,1);
  fitFun->SetParLimits(2,1,1);
  fitFun->SetParLimits(3,1,1);
  fitFun->SetParLimits(4,1,1);
  
  fitFun->SetParName(0,"N1");
  fitFun->SetParName(1,"binwidth");
  fitFun->SetParName(2,"sigma ");
  fitFun->SetParName(4,"M1");
  fitFun->SetParName(3,"G1");

  fitFun->SetParName(9,"n2");
  fitFun->SetParName(10,"G2");
  fitFun->SetParName(11,"M2");

  h->Fit("ddst","QR");
  h->Fit("ddst","QR");
  fitFun->SetParLimits(0,0,1000000);
  fitFun->SetParLimits(9,0,1000000);
  fitFun->SetParLimits(4,0.4,0.5);
  fitFun->SetParLimits(11,0.4,0.5);
  h->Fit("ddst","QR");
  h->Fit("ddst","QR");
  h->Fit("ddst","EMVR");
  h->Fit("ddst","EMVR");
  delete fitFun;
}
