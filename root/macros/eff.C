#include "TH1F.h"
#include "TMath.h"

//
//
// this program fills histogram with the gaussians
// parametes histigram pointer, mass, sigma, number of events 
//
//


void exp(TH1F* h, double m0=1.862, double pmax=0.006,int n=922) { 
	int i=0;
	while (i<n) { 
		double prb    = -log(rand()*pmax)*m0;
		if (prb<=gauss) { 
			h->Fill(m0+mass);
			i++;
		}
	}
}
