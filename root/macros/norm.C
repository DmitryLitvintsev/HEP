#include "stdio.h"
#include "TH1F.h"

void norm(TH1F* h, double scale=1) { 
	double sum = 0;
	for(int i=0;i<h->GetNbinsX();i++) {
		sum +=  h->GetBinContent(i+1);
	}

	sum /= scale;
	printf("%8.6e\n",sum);

	for(int i=0;i<h->GetNbinsX();i++) {
		double d  = h->GetBinContent(i+1);
		double ed = h->GetBinError(i+1);

		d /= sum;
		ed /= sum;
		h->SetBinContent(i+1,d);
		h->SetBinError(i+1,ed);
	}
}
