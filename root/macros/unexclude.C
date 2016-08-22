
#include "TH1.h"

void unexclude(TH1F* h, double m_low=0.37, double m_high=0.5) { 
	for ( int i=0; i<=h->GetNbinsX();i++) { 
		if ( h->GetBinCenter(i+1)>m_low&& 
		     fabs(h->GetBinCenter(i+1)<m_high)) { 
			h->SetBinError(i+1,sqrt(h->GetBinContent(i+1)));
		}
	}
}
