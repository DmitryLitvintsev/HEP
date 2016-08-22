#include "TMath.h"
#include <iostream>


using namespace std;

long  fact(long n) { 
	if ( n == 0 ) return 1;
	if ( n == 1 ) return 1;
	return fact(n-1)*n;
}

double Binomial(
	const int& s,
	const int& N, 
	const double& prb) {
	double p=0.;
	double q=1.-prb;
	for (int i=0;i<s;i++) { 
		p+=TMath::Binomial(N,i)*pow(q,N-i)*pow(prb,i);
		cout << p << endl;
	}
	return 1.-p;
}

	
