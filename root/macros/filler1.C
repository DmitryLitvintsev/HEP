#include "TH1F.h"
#include "TMath.h"
#include "TList.h"
#include "TCut.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include "TVector3.h"


//
//
// this program fills histogram with the gaussians
// parametes histigram pointer, mass, sigma, number of events 
//
//


void filler1(TH1F* h, double m0=1.321, double s1=0.006, double s2=0.0018, int n=922) { 


	int i=0;
	while (i<n) { 
		double mass   = (1.-2.*(double(rand())/RAND_MAX))*s1*10.;
		
		double r     = 0.647861;
		double nsig1 = mass / s1; 
		double nsig2 = mass / s2; 
		double prb    = (double(rand())/RAND_MAX);
		
		double gauss = r/(1.+r)*exp(-0.5*nsig1*nsig1)
			+1./(1.+r)*exp(-0.5*nsig2*nsig2);
		

		if (prb<=gauss) { 
			h->Fill(m0+mass);
			i++;
		}
	}
	double alpha =  r/(1.+r);
	double beta  =  1./(1.+r);
	double ss = alpha*s1 + beta*s2;
	printf("%8.6e\n",ss);
	double ss = sqrt(alpha*s1*s1 + beta*s2*s2);
	printf("%8.6e\n",ss);
	double ss = sqrt(alpha*alpha*s1*s1 + beta*beta*s2*s2);
	printf("%8.6e\n",ss);
}


