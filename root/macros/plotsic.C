#include "TH1F.h"
#include "TFrame.h"
#include "TLatex.h"
#include "TPad.h"
#include "decorate.C"

int plotsic(TH1*h, TH1* h1) {
	h->SetMaximum(25);
	h->SetAxisRange(0.136,0.22);
	h1->SetAxisRange(0.136,0.22);
	const char txt[] = h->GetTitle();
	h->SetTitle("");


	TFrame* frame = (TFrame*) gPad->GetListOfPrimitives()->FindObject("TFrame");

	double padWidth  = gPad->GetWNDC();
	double padHeight = gPad->GetHNDC();
	double xlow      = gPad->GetXlowNDC();
	double ylow      = gPad->GetYlowNDC();

	double sX     = 0.6;
	double sY     = 0.75;
	double width  = 0.07*padWidth;
	double height = 0.07*padHeight;

	double x1   = xlow+padWidth*sX;
	double x2   = x1 + width;
	double y11  = ylow+padHeight*sY;
	double y12  = y11+height;
	
	double y21 = y12 + 0.2*height;
	double y22 = y21 + height;

	TPave* p1 = new TPave(x1,y11,x2,y12,1,"NDC");
	TPave* p2 = new TPave(x1,y21,x2,y22,1,"NDC");


	
	p1->SetBorderSize(1);
	p1->SetFillColor(h1->GetFillColor());
	p1->SetLineColor(h1->GetLineColor());

	p2->SetLineWidth(2);
	p2->SetBorderSize(1);
	p2->SetFillColor(0);
	p2->SetLineColor(1);
	
	
	TLatex* txt1 = new TLatex(x2+0.01,y12-0.5*height,"#Lambda_{c} sideband");
	txt1->SetTextFont(42);
	TLatex* txt2 = new TLatex(x2+0.01,y22-0.5*height,"#Lambda_{c} sideband");
	txt2->SetTextFont(42);
	txt1->SetTextSize(0.01);
	
	txt1->SetNDC();
	txt2->SetNDC();
	
	txt1->SetTextAlign(12);
	txt2->SetTextAlign(12);
	txt1->SetTextSize(0.065);
	txt2->SetTextSize(0.065);
	
	h->Draw();
	decorate(h,"[GeV/c^{2}]","N / 1.5 MeV/c^{2}","M(#Lambda_{c}#pi^{-})-M(#Lambda_{c})");
	p1->Draw();
//	p2->Draw();
//	txt1->Draw();
	txt1->Draw();
	h1->Draw("same");
        gPad->RedrawAxis();
	h->Draw("same");
	return 0;
	
}
