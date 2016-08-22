#include "TH1.h"

int divide(TH1* one, TH1* two, TH1* result) { 
	for (int i=1;i<=one->GetNbinsX(); i++) { 
		if ( two->GetBinContent(i)<=0) continue;
		if ( one->GetBinContent(i)<=0) continue;
		double d1  = one->GetBinError(i);
		double d2  = two->GetBinError(i);

		double a1  = one->GetBinContent(i);
		double a2  = two->GetBinContent(i);

		double e  = a1 / a2;
		double de = e*sqrt(d1*d1/a1/a1+d2*d2/a2/a2);

		result->SetBinContent(i,e);
		result->SetBinError(i,de);
		
	}
	return 0;
}
