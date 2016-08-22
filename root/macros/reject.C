
#include "TH1.h"

void reject(TH1F* h, double mass=1.5319, double sigma=0.01) { 
	for ( int i=0; i<=h->GetNbinsX();i++) { 
		if ( fabs(h->GetBinCenter(i+1)-mass)<3.*sigma) { 
			h->SetBinError(i+1,0);
		}
	}

}
 


