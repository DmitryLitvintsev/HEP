#include "TH1F.h"
#include "TMath.h"
#include "TF1.h"
#include <iostream>
#include "/cdf/atom/home/litvinse/root_macros/filler.C"
#include "/cdf/atom/home/litvinse/root_macros/fit.C"
#include "/cdf/atom/home/litvinse/root_macros/penta.C"

using namespace std;

void ul(TH1F* hh) { 
	hh->Sumw2();
	TH1F* tmp =  new TH1F("tmp",hh->GetTitle(),
			      hh->GetNbinsX(),
			      hh->GetXaxis()->GetXmin(),
			      hh->GetXaxis()->GetXmax());	

	penta(hh,4,1.4,2.4,1.862,0.008,0);

	double chi2_0 = hh->GetFunction("pentafit")->GetChisquare();
	int  NDF_0    = hh->GetFunction("pentafit")->GetNDF();
	
	double prb_0 = TMath::Prob(chi2_0,NDF_0);
	double prb = prb_0;
	double r = prb / prb_0;
	cout << chi2_0 << " " << NDF_0 << endl;
	int n = 0;
	while ( r > 0.1 ) { 
		filler(tmp,1.862,0.007,n++);
		tmp->Sumw2();
		TH1F h100 = (*tmp)+(*hh);
		penta(&h100,4,1.4,2.4,1.862,0.008,0);
		double chi2_n = h100.GetFunction("pentafit")->GetChisquare();
		int  NDF_n    = h100.GetFunction("pentafit")->GetNDF()-1;
		double prb_n = TMath::Prob(chi2_n,NDF_0-1);
		r = prb_n / prb_0;
		tmp->Reset();
	}
	cout << r << " " << n << endl;
		 
}
