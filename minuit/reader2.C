#include "TF1.h"
#include "TH1F.h"
#include "TList.h"
#include "TCut.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include "TVector3.h"
#include "TPolyLine.h"

#include "/cdf/atom/home/litvinse/root_macros/decorate.C"
#include "TStyle.h"
#include "TPad.h"
#include "TLatex.h"
#include "TArrow.h"
#include "TPolyLine.h"
#include "TLine.h"
#include "utils.h"
#include "fit.C"

#include <algorithm>

//----------------------------------------------------------------------
//
// this one plots UL vs mass D0 p 
//
//
//----------------------------------------------------------------------



 
void reader2(TH1F* h, TF1* f, char* fname) 
{
	
	f = new TF1("xist",xist,1.4,2.2,10);

	FILE* fp = fopen(fname,"r");
	float par,eparab,eminus,eplus;
	char pname[10];
	int npar,ii;
	while(fscanf(fp,"%d %s %f %f %f %f %d",&npar,&pname[0],&par,&eparab,&eminus,&eplus,&ii)!=EOF) { 
		printf("%d %s %f %f %f %f\n",npar,pname,par,eparab,eminus,eplus);
		f->SetParameter(npar,(double)par);	
		f->SetParError(npar,(double)eparab);
		f->SetParName(npar,pname);
	}
	
	f->SetParameter(0,h->GetBinWidth(1));
	f->SetParameter(0,h->GetBinWidth(1));
	f->SetParLimits(0,1,1);


	double chi2 = 0;
        int ndof = 0;
        for(int i=0; i<h->GetNbinsX();i++) { 
                double x      = h->GetBinCenter(i+1);
                double y      = h->GetBinContent(i+1);
		if(y<0.5) continue;
		if (x>2.2) continue;
                double dy     = sqrt(y);
                double theory = f->Eval(x);
                chi2 += ( y - theory ) * ( y - theory ) / dy / dy;
                ndof ++;
        }
	ndof -= npar;

        double chi2_per_ndf = chi2 / (double)(ndof);
        double prob         = TMath::Prob(chi2,ndof);

        cout << chi2 << " " << ndof << endl;

        char txt[50];
        sprintf(txt,"#chi^{2}/ndf=%2.2f/%d=%2.2f Prob=%2.2f",chi2,ndof,chi2_per_ndf,prob);

//        printf("chi2/ndf=%d/%d=%2.2f Prob =%2.2f",chi2,ndof,chi2_per_ndf,prob);

        gStyle->SetOptFit(0);
        gStyle->SetOptStat(0);
        h->GetXaxis()->SetNdivisions(510);
        h->GetYaxis()->SetNdivisions(505);
        h->SetTitle("CDF Run II Preliminary");
	h->SetAxisRange(1.4,2.2);
        h->Draw();
        double padWidth  = gPad->GetWNDC();
        double padHeight = gPad->GetHNDC();
        double xlow      = gPad->GetXlowNDC()+0.2*padWidth;
        double ylow      = gPad->GetYlowNDC();
        
        double sX     = 0.2;
        double sY     = 0.8;

        double x = xlow+sX*padWidth;
        double y = ylow+sY*padHeight;
        TLatex* ltxt = new TLatex(x,y,txt);     
        ltxt->SetNDC();
        ltxt->SetTextAlign(11);
        ltxt->SetTextSize(0.05);
        decorate(h,"[GeV/c^{2}]","N / 5 MeV/c^{2}"," M(#Xi #pi^{+})");
        f->SetLineColor(4);
        f->Draw("same");
        ltxt->Draw("same");



	
}

