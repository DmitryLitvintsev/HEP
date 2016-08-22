#include "TH1.h"
#include "TMath.h"
#include <iostream>

using namespace std;

void compare(TH1* h1, TH1* h2) { 
	double chi2 = 0.;
	int      np = 0;
	for (int i=0;i<h1->GetNbinsX();i++) { 
		double a1  = h1->GetBinContent(i+1);
		double da1 = h1->GetBinError(i+1);
		double a2  = h2->GetBinContent(i+1);
		double da2 = h2->GetBinError(i+1);
		if (a1<=0||a2<=0) continue;
		chi2 += (a1-a2)*(a1-a2)/(da1*da1+da2*da2);
		np++;
	}
	printf("%8.6e / %d %8.6e %8.6e\n",chi2, np, chi2/(double)np,TMath::Prob(chi2,np));
	
}
