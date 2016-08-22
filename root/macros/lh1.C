#include "TH1F.h"
#include "TList.h"
#include "TCut.h"
#include "TMath.h"
#include <iostream>
#include <vector>
#include "/cdf/atom/home/litvinse/root_macros/fit.C"
#include "/cdf/atom/home/litvinse/root_macros/penta.C"
#include <fstream>
#include <iterator>

using namespace std;

void lh1(TH1F* hh, TH1F* h2,int n=100) {
	vector<double> chi2s;
	cout << n << endl;
	for (int i=0; i<n; i++) { 
		penta(hh,5,1.6,2.4,1.862,0.008,i);
//		penta(hh,5,-1,2.4,1.862,0.008,i);
		chi2s.push_back(hh->GetFunction("pentafit")->GetChisquare());
	}
	double minchi2 = 1.e+10;
	for(vector<double>::const_iterator i=chi2s.begin();
	    i!=chi2s.end();i++) { 
		if ((*i)<minchi2) minchi2=(*i);
	}

	vector<double> lhs;
	double sum=0;

	for(vector<double>::const_iterator i=chi2s.begin();
	    i!=chi2s.end();i++) { 
		double chi2 = (*i)-minchi2;
		double lh   = exp(-0.5*chi2);
		sum += lh;
		lhs.push_back(lh);
	}
	
	int j = 1;
	double r = 0;
	for(vector<double>::const_iterator i=lhs.begin();
	    i!=lhs.end();i++) { 
		h2->SetBinContent(j,(*i));
		r += (*i)/sum;
		if (fabs(r-0.93)<0.05) {
			cout << r << " " << j-1 << endl;
		}
		j++;
	}
}
