#include "TMath.h"
#include <iostream>


using namespace std;

long  fact(long n) { 
	if ( n == 0 ) return 1;
	if ( n == 1 ) return 1;
	return fact(n-1)*n;
}

double prob(const double& B, 
	    const int& N) {
	double p=0.;
	for (int i=0;i<N-1;i++) { 
		p+=pow(B,i)/TMath::Factorial(i)*exp(-B);
	}
	return 1.-p;
}

	
