#include "TH1F.h"
#include <iostream>

using namespace std;

void integral4(TH1F* hh) {
	double sum = 0;
	for (int i=0; i<hh->GetNbinsX();i++) { 
		if (hh->GetBinLowEdge(i+1)<0) continue;
		sum += hh->GetBinContent(i+1);
	}

	double s = 0;
	for (int i=0; i<hh->GetNbinsX();i++) { 
		if (hh->GetBinLowEdge(i+1)<0) continue;
		s += hh->GetBinContent(i+1);
		double r = s / sum;
		if ( fabs(r-0.9)<0.01 ) cout << r << " " <<  hh->GetBinLowEdge(i+1) << endl;
		if ( fabs(r-0.95)<0.01 ) cout << r << " " <<  hh->GetBinLowEdge(i+1) << endl;
	}
	
}
