#include "TH1F.h"
#include <iostream>

using namespace std;

double get_integral(TH1F* hh) {
	double sum = 0;

	for (int i=0; i<hh->GetNbinsX();i++) { 
		sum += hh->GetBinContent(i+1);
	}

	return sum;
	
}
