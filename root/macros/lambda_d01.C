#include "TH1F.h"
#include "TList.h"
#include "TCut.h"
#include "TF1.h"
#include "TMinuit.h"
#include "hyperon.h"
#include "hyperon.C"
#include "/cdf/atom/home/litvinse/root_macros/fit.C"
#include "/cdf/atom/home/litvinse/root_macros/fitter2.C"
#include <iostream>

using namespace std; 

//
// Efficiency vs d0 cuts on lambda tracks, SVX xi
//
//

void lambda_d01(TList* right,TList* wrong, TList* result,hyperon& xi)
{
        Int_t nbins = 40;
        Double_t var_max = 0.5;
        Double_t var_min = 0.;

	TH1F* h600 = new TH1F("significance_vs_d0","significance_vs_d0",nbins,var_min,var_max);
	double delta = (var_max-var_min)/(double)nbins;

       for (int i=0; i<nbins; i++) { 
		char txt[10];
		sprintf(txt,"rs%d",i);
		TH1F* ptr1 = new TH1F(txt,txt,80,1.24,1.4);
		sprintf(txt,"ws%d",i);
		TH1F* ptr2 = new TH1F(txt,txt,80,1.24,1.4);

		right->Add(ptr1);
		wrong->Add(ptr2);
	}

	for (int j=0; j<xi.fChain->GetEntries(); j++) {
		if (j%100000==0) cout << " doing " << j << endl;
		xi.GetEntry(j);
		if (xi._1_v0_vrt<xi._1_h_vrt+1)    continue;
		if (xi._1_nsvx<0.5)                continue;						   

		for (int i=0; i<nbins; i++) { 
// 			if (fabs(xi._1_kink_d0)>delta*(double)(i)&&fabs(xi._1_pi_d0)>delta*(double)(i))
			if (fabs(xi._1_pr_d0)>delta*(double)(i) &&fabs(xi._1_pi_d0)>delta*(double)(i)) { 
				if (xi._1_q1*xi._1_q2<0) {
					((TH1F*)right->At(i))->Fill(xi._1_h_mass);
				}
				else { 
					((TH1F*)wrong->At(i))->Fill(xi._1_h_mass);	
				}
			}
		}
	}
	 
	for (int i=0; i<nbins; i++) { 
		TH1F* ptr1 = (TH1F*)right->At(i);
		TH1F* ptr2 = (TH1F*)wrong->At(i);
		
		char txt[10];
		sprintf(txt,"diff_%d",i);
		
		TH1F* sum =  new TH1F(txt,txt,80,1.24,1.4);

		sum->SetAxisRange(1.29,1.35);
		ptr1->Sumw2();
		ptr2->Sumw2();
		sum->Add(ptr1,ptr2,1.,-1.);

		fitter2(sum,3,1.29,1.35);
		cout << "fitting " << sum->GetTitle() << " " << gMinuit->fCstatu  << endl;

		double signal = sum->GetFunction("gp2")->GetParameter(0);
		double sigma  = sum->GetFunction("gp2")->GetParError(0);

		double r = signal / sigma;
		result->Add(sum);


		h600->SetBinContent(i+1,r);
		h600->SetBinError(i+1,0);
	}
	result->Add(h600);
	       
}

