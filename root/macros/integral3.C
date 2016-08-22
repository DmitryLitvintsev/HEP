#include "TH1F.h"
#include <iostream>

using namespace std;

void integral3(TH1F* hh) {
	double sum = 0;
	for (int i=0; i<hh->GetNbinsX();i++) { 
		if (hh->GetBinLowEdge(i+1)<0) continue;
		sum += hh->GetBinContent(i+1);
	}
	cout << sum << endl;

	double s = 0;
	for (int i=0; i<hh->GetNbinsX();i++) { 
		if (hh->GetBinLowEdge(i+1)<0) continue;
		s += hh->GetBinContent(i);
		double r = s / sum;
		cout << r << " " << i << " " <<  hh->GetBinLowEdge(i+1) << endl;
	}
	
}
