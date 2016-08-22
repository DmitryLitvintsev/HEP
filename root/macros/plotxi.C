#include "TH1F.h"
#include "TFrame.h"
#include "TLatex.h"
#include "TPad.h"

int plotlamcst(TH1*h) {
	TFrame* frame = (TFrame*) gPad->GetListOfPrimitives()->FindObject("TFrame");

	double padWidth  = gPad->GetWNDC();
	double padHeight = gPad->GetHNDC();
	double xlow      = gPad->GetXlowNDC();
	double ylow      = gPad->GetYlowNDC();

	double sX     = 0.6;
	double sY     = 0.7;
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
	p1->SetFillColor(1);
	p1->SetFillStyle(3001);
	p1->SetLineWidth(2);
	p1->SetLineColor(1);

	p2->SetLineWidth(2);
	p2->SetBorderSize(1);
	p2->SetFillColor(0);
	p2->SetLineColor(1);
	
	
	TLatex* txt1 = new TLatex(x2+0.01,y22-0.5*height,"#Lambda#pi^{-}");
	TLatex* txt2 = new TLatex(x2+0.01,y12-0.5*height,"#Lambda#pi^{+}");
	
	txt1->SetNDC();
	txt2->SetNDC();
	
	txt1->SetTextAlign(12);
	txt2->SetTextAlign(12);
	txt1->SetTextSize(0.065);
	txt2->SetTextSize(0.065);
	
	p1->Draw();
	p2->Draw();
	txt1->Draw();
	txt2->Draw();
	return 0;
	
}
