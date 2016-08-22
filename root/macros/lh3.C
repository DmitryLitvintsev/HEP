#include "TH1F.h"
#include "TList.h"
#include "TCut.h"
#include "TMath.h"
#include "TAxis.h"
#include <iostream>
#include <vector>
#include "/cdf/atom/home/litvinse/root_macros/fit.C"
#include "/cdf/atom/home/litvinse/root_macros/penta3.C"
#include <fstream>
#include <iterator>

using namespace std;

//
// this macro was used to produce LH scan for (Xi-,pi+) BGW + SQRT(x)*Pol
//

void lh3(TH1F* hh, TH1F* h2) {
	TAxis* a  = h2->GetXaxis();
	
	vector<double> chi2s;
	for (int i=1; i<=h2->GetNbinsX(); i++) { 
		cout << " doing " << i << endl;
		int nevents = h2->GetBinLowEdge(i);
		penta3(hh,4,1.4,2.4,1.862,0.008,nevents);
		chi2s.push_back(hh->GetFunction("pentafit1")->GetChisquare());
	}
	double minchi2 = 1.e+10;

	for(vector<double>::const_iterator it=chi2s.begin();
	    it!=chi2s.end();it++) { 
		if ((*it)<minchi2) minchi2=(*it);
	}

	vector<double> lhs;
	double sum=0;

	for(vector<double>::const_iterator it=chi2s.begin();
	    it!=chi2s.end();it++) { 
		double chi2 = (*it)-minchi2;
		double lh   = exp(-0.5*chi2);
		sum += lh;
		lhs.push_back(lh);
	}
	
	int j = 1;
	double r = 0;
	for(vector<double>::const_iterator it=lhs.begin();
	    it!=lhs.end();it++) { 
		h2->SetBinContent(j,(*it));
		r += (*it)/sum;
		if (fabs(r-0.93)<0.1) {
			cout << r << " " << j << endl;
		}
		j++;
	}
}
