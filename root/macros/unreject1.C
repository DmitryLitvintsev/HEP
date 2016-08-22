#include "TH1.h"

void unreject1(TH1F* h, double low=0.37, double high=0.5) {
	for ( int i=0; i<=h->GetNbinsX();i++) { 
		if (h->GetBinCenter(i+1)<high&&h->GetBinCenter(i+1)>low) { 
			h->SetBinError(i+1,sqrt(h->GetBinContent(i+1)));
		}
	}
	
}
 
