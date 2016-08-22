#include "TH1F.h"
#include "TLatex.h"
#include "TPad.h"
#include "TFrame.h"
#include "TF1.h"


int fitter1decorate(TH1*h, 
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

	double x = xlow+sX*padWidth;
	double y = ylow+sY*padHeight;

	TLatex* txt = new TLatex(x,y,title);
	txt->SetNDC();
	txt->SetTextSize(0.08);
	txt->SetTextAlign(21);
	txt->SetTextFont(42);
	txt->Draw();

	sX     = 0.6;
	sY     = 0.8;

        double height = 0.06*padHeight;

        double x1   = xlow+padWidth*sX;
        double y1   = ylow+padHeight*sY;

	TF1* f = (TF1*)h->GetFunction("gp");
	

	char text[50];
	sprintf(text,"Yield=%d#pm%d",(int)f->GetParameter(0)+0.5,
		(int)f->GetParError(0)+0.5);
	
	TLatex* t1 = new TLatex(x1,y1,text);
	t1->SetNDC();
	t1->SetTextSize(0.05);
	t1->SetTextFont(42);
	t1->Draw();

        y1  -= height;

	sprintf(text,"M=(%4.1f#pm%4.1f)MeV/c^{2}",
		f->GetParameter(3)*1000.,
		f->GetParError(3)*1000.);
	
	TLatex* t2 = new TLatex(x1,y1,text);
	t2->SetNDC();
	t2->SetTextSize(0.05);
	t2->SetTextFont(42);
	t2->Draw();

        y1  -= height;

	sprintf(text,"#sigma=(%4.1f#pm%4.1f)MeV/c^{2}",
		f->GetParameter(2)*1000.,
		f->GetParError(2)*1000.);
	
	TLatex* t3 = new TLatex(x1,y1,text);
	t3->SetNDC();
	t3->SetTextSize(0.05);
	t3->SetTextFont(42);
	t3->Draw();


        y1  -= height;

	sprintf(text,"prob=%3.1f%\n",f->GetProb()*100.);
	TLatex* t4 = new TLatex(x1,y1,text);
	t4->SetNDC();
	t4->SetTextSize(0.05);
	t4->SetTextFont(42);
	t4->Draw();

        x1   = xlow+padWidth*0.8;
	y1   = ylow+padHeight*0.92;

	TLatex* t5 = new TLatex(x1,y1,extra);
	t5->SetNDC();
	t5->SetTextSize(0.06);
	t5->SetTextFont(42);
	t5->Draw();
	return 0;
	
}
