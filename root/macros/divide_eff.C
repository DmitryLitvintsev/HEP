#include "TH1.h"

int divide_eff(TH1* one, TH1* two, TH1* result) { 
	//
	// one / two = result
	//
	for (int i=1;i<=one->GetNbinsX(); i++) { 
		if ( two->GetBinContent(i)<=0) continue;
		if ( one->GetBinContent(i)<=0) continue;
		double e  = one->GetBinContent(i) / two->GetBinContent(i);
		double de = sqrt((1.-e)*e / one->GetBinContent(i));
		result->SetBinContent(i,e);
		result->SetBinError(i,de);
	}
	return 0;
}
