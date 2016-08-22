#include "TH1F.h"
#include "TLatex.h"
#include "TPad.h"


int decorate(TH1*h, 
	     const char* xTitle,
	     const char* yTitle,
	     const char* title,
	     const char* extra="") { 
	h->GetXaxis()->SetTitle(xTitle);
	h->GetYaxis()->SetTitle(yTitle);
	h->GetXaxis()->SetTitleSize(0.06);
	h->GetYaxis()->SetTitleSize(0.06);
	h->GetYaxis()->SetTitleOffset(1.2);
	h->GetXaxis()->SetTitleOffset(1.0);
	h->GetXaxis()->SetLabelSize(0.05);
	h->GetYaxis()->SetLabelSize(0.05);
	h->GetXaxis()->SetLabelOffset(0.007);
	h->GetYaxis()->SetLabelOffset(0.007);

	h->GetYaxis()->SetLabelFont(42);
	h->GetYaxis()->SetTitleFont(42);

	h->GetXaxis()->SetLabelFont(42);
	h->GetXaxis()->SetTitleFont(42);

	double padWidth  = gPad->GetWNDC();
	double padHeight = gPad->GetHNDC();
	double xlow      = gPad->GetXlowNDC();
	double ylow      = gPad->GetYlowNDC();

	
	double sX     = 0.55;
	double sY     = 0.02;
//	double sY     = 0.15;

	double x = xlow+sX*padWidth;
	double y = ylow+sY*padHeight;

	TLatex* txt = new TLatex(x,y,title);
	txt->SetNDC();
	txt->SetTextSize(0.08);
	txt->SetTextAlign(21);
	txt->SetTextFont(42);
	txt->Draw();
        double x1   = xlow+padWidth*0.8;
	double y1   = ylow+padHeight*0.92;

	TLatex* t5 = new TLatex(x1,y1,extra);
	t5->SetNDC();
	t5->SetTextSize(0.06);
	t5->SetTextFont(42);
	t5->Draw();
	return 0;
	
}
