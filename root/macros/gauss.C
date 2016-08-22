#include <math>
#include <stdio.h>
#include <stdlib.h>


double gauss(const double& mean=0, const double& sigma=1) {
	while (1) { 
		double mass   = (1.-2.*(double(rand())/RAND_MAX))*sigma*0.5;
		double prb    = (double(rand())/RAND_MAX);
		double gauss  =  exp(-0.5*mass*mass/sigma/sigma);
		if (prb<=gauss) { 
			return mass+mean;
		}
	}
}
