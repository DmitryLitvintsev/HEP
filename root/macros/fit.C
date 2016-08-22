#include "TMath.h"
Double_t Gauss(Double_t *x, Double_t mass, Double_t sigma) {
	double arg = (x[0]-mass)/sigma;
	return exp(-0.5*arg*arg)/sqrt(2.*TMath::Pi())/sigma;
}

Double_t bu_fit(Double_t *x, Double_t *par) {
	static const double mk2 = 2.43716980329e-1;
	static const double mpi2  = 1.94798351452323965e-2;
	double alpha = x[1];

//	if (alpha<0.5) return 0;

	double mm = par[3]*par[3] + (1.+alpha)*(mk2-mpi2);
	mm=sqrt(mm);

	double fk[3] = { 0. , 0.149 , 0.774 };
	fk[0] = 1.-fk[1]-fk[2];

	double fpi[3] = { 0. , 0.134 , 0.772 };
	fpi[0] = 1.-fpi[1]-fpi[2];

	double fbg[3] = { 0. , 0.130 , 0.846 };
	fbg[0] = 1.-fbg[1]-fbg[2];

	double lambda_k[] =  { 0.600, 4.629, 1.275};
	double lambda_pi[] = { 0.631, 5.013, 1.230}; 
	double lambda_bg[] = { 0.726, 6.540, 1.978};

	double ak = 0.527;
	double api = 0.522;
	double abg = 0.398;

	double hk=0.;
	double hpi=0.;
	double hbg=0.;
	for (int i=0;i<3;i++) { 
		hk+=fk[i]*(alpha-ak)*exp(-lambda_k[i]*alpha);
		hpi+=fpi[i]*(alpha-api)*exp(-lambda_pi[i]*alpha);
		hbg+=fbg[i]*(alpha-abg)*(alpha-abg)*(alpha-abg)*
			exp(-lambda_bg[i]*alpha);
	}
		
	
	return par[0]*(par[1]*Gauss(x,par[3],par[2])*hk+
		       par[4]*Gauss(x,mm,par[2])*hpi+
		       (par[5]+par[6]*x[0])*hbg);
}
		       

Double_t bu_fit2(Double_t *x, Double_t *par) {
	static const double mk2 = 2.43716980329e-1;
	static const double mpi2  = 1.94798351452323965e-2;
	double alpha = x[1];

	if (alpha<0.5) return 0;

	double n1   = par[1];
	double s1   = par[2];
	double s2   = par[3];
	double r    = par[4];
	double mass = par[5];
	double n2   = par[6];
	double a0   = par[7];
	double a1   = par[8];
	double a2   = par[9];

 	double mm = mass*mass + (1.+alpha)*(mk2-mpi2);
 	mm=sqrt(mm);

 	double fk[3] = { 0. , 0.149 , 0.774 };
 	fk[0] = 1.-fk[1]-fk[2];

 	double fpi[3] = { 0. , 0.134 , 0.772 };
 	fpi[0] = 1.-fpi[1]-fpi[2];

 	double fbg[3] = { 0. , 0.130 , 0.846 };
 	fbg[0] = 1.-fbg[1]-fbg[2];

 	double lambda_k[] =  { 0.600, 4.629, 1.275};
 	double lambda_pi[] = { 0.631, 5.013, 1.230}; 
 	double lambda_bg[] = { 0.726, 6.540, 1.978};

 	double ak = 0.527;
 	double api = 0.522;
 	double abg = 0.398;

 	double hk=0.;
 	double hpi=0.;
 	double hbg=0.;
 	for (int i=0;i<3;i++) { 
 		hk+=fk[i]*(alpha-ak)*exp(-lambda_k[i]*alpha);
 		hpi+=fpi[i]*(alpha-api)*exp(-lambda_pi[i]*alpha);
 		hbg+=fbg[i]*(alpha-abg)*(alpha-abg)*(alpha-abg)*
 			exp(-lambda_bg[i]*alpha);
 	}
		
	
 	return par[0]*(n1*
 		       (r/(1.+r)*Gauss(x,mass,s1)+Gauss(x,mass,s2)/(1.+r))*hk+
 		       n2*Gauss(x,mm,s2)*hpi+
 		       (a0+a1*x[0]+a2*x[0]*x[0])*hbg);
}
		       

Double_t eff(Double_t *x, Double_t *par) {
        return par[0]/(1.+exp(par[1]*x[0]));
}

Double_t Threshold(Double_t *x, Double_t *par) {
	return par[0]*(1.-TMath::Erf((x[0]-par[1])/par[2]));
}
Double_t GauExp(Double_t *x, Double_t *par) {
	//
	// Convolution of Gaussian and Exponential
	//
	double n       = par[0];
	double sigma  = fabs(par[1]);
	double m0     = par[2];
	double lambda = par[3];
	double xx     = x[0];
	return n/lambda*exp(pow(sigma/lambda,2)-(xx-m0)/lambda)*(1.-TMath::Freq(sigma/lambda-(xx-m0)/sigma));
}

Double_t r01(Double_t *x, Double_t *par) { 
 	double n      = 4.89374e-01;
 	double sigma  = 4.39990e-01;
 	double m0     = 2.82127e+00;
 	double lambda = 1.80454e+00;
	double xx     = x[0];
	double f_data =  n/lambda*exp(pow(sigma/lambda,2)-(xx-m0)/lambda)*(1.-TMath::Freq(sigma/lambda-(xx-m0)/sigma));
	n = 4.64590249006894052e-01;
	sigma = 6.38046697785285111e-01;
	m0 = 2.98945720305020579e+00;
	lambda = 1.84178303403883681e+00;

	double f_mc =  n/lambda*exp(pow(sigma/lambda,2)-(xx-m0)/lambda)*(1.-TMath::Freq(sigma/lambda-(xx-m0)/sigma));
	return f_data / f_mc;	
	
	
}
Double_t r1(Double_t *x, Double_t *par)  {

//
// ratio of data to MC 
//
	
// data:
 	double n      = 4.89374e-01;
 	double sigma  = 4.39990e-01;
 	double m0     = 2.82127e+00;
 	double lambda = 1.80454e+00;
	
	
	
	double xx     = x[0];
	double f_data =  n/lambda*exp(pow(sigma/lambda,2)-(xx-m0)/lambda)*(1.-TMath::Freq(sigma/lambda-(xx-m0)/sigma));
	
// 1530:

	n      = 4.63205e-01;
	sigma  = 6.37727e-01;
	m0     = 2.96660e+00;
	lambda = 1.83925e+00;

// 	n      = 4.71900983772074789e-01;
// 	sigma  = 6.41132497037125160e-01;
// 	m0     = 2.96576953646857522e+00;
// 	lambda = 2.14219176741575978e+00;

	double f_mc =  n/lambda*exp(pow(sigma/lambda,2)-(xx-m0)/lambda)*(1.-TMath::Freq(sigma/lambda-(xx-m0)/sigma));
	return f_data / f_mc;	
}


Double_t r10(Double_t *x, Double_t *par)  {
// 
// fit to pythia w/ no cuts
// 
	double n      = 4.89374e-01;
	double sigma  = 4.39990e-01;
	double m0     = 2.82127e+00;
	double lambda = 1.80454e+00;
	double xx     = x[0];
	double f_data =  n/lambda*exp(pow(sigma/lambda,2)-(xx-m0)/lambda)*(1.-TMath::Freq(sigma/lambda-(xx-m0)/sigma));
	
// 1530:
	n      = 4.64590e-01;
	sigma  = 6.38047e-01;
	m0     = 2.98946e+00;
	lambda = 1.84178e+00;
	double f_mc =  n/lambda*exp(pow(sigma/lambda,2)-(xx-m0)/lambda)*(1.-TMath::Freq(sigma/lambda-(xx-m0)/sigma));
	return f_data / f_mc;	
}


Double_t r2(Double_t *x, Double_t *par)  {
//
// make flat 1530 look like data 
//
	double n      = 4.89374e-01;
	double sigma  = 4.39990e-01;
	double m0     = 2.82127e+00;
	double lambda = 1.80454e+00;
	double xx     = x[0];
	double f_data =  n/lambda*exp(pow(sigma/lambda,2)-(xx-m0)/lambda)*(1.-TMath::Freq(sigma/lambda-(xx-m0)/sigma));
	
// 1530:
	n      = 4.71351e-01;
	sigma  = 5.66420e-01;
	m0     = 3.27838e+00;
	lambda = 1.78820e+00;
	double f_mc =  n/lambda*exp(pow(sigma/lambda,2)-(xx-m0)/lambda)*(1.-TMath::Freq(sigma/lambda-(xx-m0)/sigma));
	return f_data / f_mc;	
}

Double_t r3(Double_t *x, Double_t *par)  {
//
// make flat 1530 look like pythia MC
//
	double n      = 4.63205e-01;
	double sigma  = 6.37727e-01;
	double m0     = 2.96660e+00;
	double lambda = 1.83925e+00;
	double xx     = x[0];
	double f_data =  n/lambda*exp(pow(sigma/lambda,2)-(xx-m0)/lambda)*(1.-TMath::Freq(sigma/lambda-(xx-m0)/sigma));
	
// 1530:
	n      = 4.53338e-01;
	sigma  = 6.92564e-01;
	m0     = 3.53146e+00;
	lambda = 1.69515e+00;
	double f_mc =  n/lambda*exp(pow(sigma/lambda,2)-(xx-m0)/lambda)*(1.-TMath::Freq(sigma/lambda-(xx-m0)/sigma));
	return f_data / f_mc;	
}

 Double_t r4(Double_t *x, Double_t *par)  {
 //
 // make 1696 look like 1530 data 
 //
 	double n      = 4.89374e-01;
 	double sigma  = 4.39990e-01;
 	double m0     = 2.82127e+00;
 	double lambda = 1.80454e+00;

 	double xx     = x[0];
 	double f_data =  n/lambda*exp(pow(sigma/lambda,2)-(xx-m0)/lambda)*(1.-TMath::Freq(sigma/lambda-(xx-m0)/sigma));
	
 // :
 	n      = 4.54424e-01;
 	sigma  = 7.00642e-01;
 	m0     = 2.91861e+00;
 	lambda = 1.83467e+00;
 	double f_mc =  n/lambda*exp(pow(sigma/lambda,2)-(xx-m0)/lambda)*(1.-TMath::Freq(sigma/lambda-(xx-m0)/sigma));
 	return f_data / f_mc;
 }


//
// BW (x) Gauss
//
double bwg(const double& x, 
	   const double& gamma, 
	   const double& sigma, 
	   const double& m0) { 
	
	double epsilon = 0.5*gamma/sigma;
	double m       = x;
	double nsigm   = (m-m0)/sigma;

	
	if (gamma<0.0001) {
		return 0.3989423/sigma*exp(-0.5*nsigm*nsigm);
	}

	int npoints  = 100;
	double left  = (-5.+nsigm)/epsilon;
	double right = (5.+nsigm)/epsilon;
	double step  = (right-left)/npoints;
	double sum   = 0.;

	for (int i=1; i<=npoints; i++) { 
		double var=left+((double)i-0.5)*step;
		sum+=exp(-0.5*(nsigm-epsilon*var)*(nsigm-epsilon*var))/(var*var+1.);
	}
	return sum*step*2./sigma/15.7496;
}

double spectrum(double *x, double *par) {
	return par[0]*pow(x[0],par[1])*exp(-par[2]*pow(x[0],par[3]));
}

double spectrum1(double *x, double *par) {
	return exp(par[0]+par[1]*x[0]+par[2]*pow(x[0],2));
}

//
// BW(x)Gauss+polynom&sqrt
//

double bgwpol(double *x, double *par) {
	double thr = .1396+2.010;
	double bin = 10.;
	double dm  = x[0]-thr;

	if ( x[0]<thr ) return 0;

	return bin*par[0]*bwg(x[0],par[2],par[1],par[3])+
		bin*par[4]*bwg(x[0],par[6],par[5],par[7])+
		sqrt(dm)*fabs(par[8])*(1.+dm*par[9]+par[10]*dm*dm);
}

double bwg0(double *x, double *par) {
	return par[0]*bwg(x[0],par[2],par[1],par[3]);
}

//!!!!!!!!!!!!!!!!!
Double_t spec1(Double_t *x, Double_t *par) {
	return par[0]*pow(x[0],par[1])*exp(par[2]*x[0]);
}

Double_t spec2(Double_t *x, Double_t *par) {
	return par[0]*Gauss(x,par[1],par[2])*exp(par[3]*pow(x[0],par[4]));
}


Double_t Polynomial(Double_t *x, Double_t* par) 
{
	double dx = x[0]-par[0];
	return par[1]+par[2]*dx+par[3]*dx*dx+par[4]*dx*dx*dx;
}

Double_t Exp(Double_t *x, Double_t* par) 
{
	double dx = x[0]-par[0];
	return exp(par[1]+par[2]*dx+par[3]*dx*dx+par[4]*dx*dx*dx);
}

Double_t Polynomial1(Double_t *x, Double_t* par) 
{
	double dx = x[0]-(2.010 + 0.93827231);
	return par[0]+par[1]*dx+par[2]*dx*dx+par[3]*dx*dx*dx;
}



Double_t lambdafit(Double_t *x, Double_t *par) {
	Double_t thr = 0.13962 + 0.93827231;
	return par[1]*(par[0]*par[3]/(1.+par[3])*Gauss(x,par[5],par[2])+
		       par[0]/(1.+par[3])*Gauss(x,par[5],par[4]))+
		sqrt(0.5*(fabs(x[0]-thr)+x[0]-thr))*
		Polynomial(x,&par[5]);
}


Double_t dst2(Double_t *x, Double_t *par) {
	Double_t thr = 2.010 + 0.93827231;
//	Double_t thr = 0.13962 ;
	return par[1]*par[0]*bwg(x[0],par[3],par[2],par[4])+
		sqrt(0.5*(fabs(x[0]-thr)+x[0]-thr))*
		Polynomial1(x,&par[5])*3.e-3;
}


Double_t lamcstfit(Double_t *x, Double_t *par) {
	if (x[0]<2.*0.1396) return 0;
	Double_t dmass = x[0]-2.*0.1396;
	return par[0]*(par[1]*bwg(x[0],par[3],par[2],par[4])+
		       par[5]*Gauss(x,par[6],par[2]))+
		par[0]*pow(dmass,par[7])*(par[8]+par[9]*dmass);
}

Double_t sigmacfit(Double_t *x, Double_t *par) {
	if (x[0]<0.1396) return 0;
	Double_t dm = x[0]-0.1396;
	return par[0]*(par[1]*bwg(x[0],0.,par[2],par[3]))+
		par[0]*pow(dm,par[4])*(par[5]+par[6]*dm+par[7]*dm*dm+par[8]*dm*dm*dm);
}





  Double_t BW(Double_t *x, Double_t mass, Double_t gamma) 
  {
	  return 0.5*gamma/TMath::Pi()/((x[0]-mass)*(x[0]-mass)+0.25*gamma*gamma);
  }


  Double_t bw(Double_t *x, Double_t* par) 
  {
	  return 0.5*par[1]/TMath::Pi()/((x[0]-par[0])*(x[0]-par[0])+0.25*par[1]*par[1]);
  }

  Double_t xifit(Double_t *x, Double_t *par) {
	  return par[0]*par[1]*Gauss(x,par[9],par[2])+
		 par[3]*par[1]*Gauss(x,par[8],par[4])+
		 par[5]*par[1]*BW(x,par[7],par[6])+
		 sqrt(0.5*(fabs(x[0]-(1.115+0.1396))+x[0]-(1.115+0.1396)))*Polynomial(x,&par[9]);
  }

  Double_t xifit0(Double_t *x, Double_t *par) {
	  return par[0]*par[1]*Gauss(x,par[3],par[2])+
		 sqrt(0.5*(fabs(x[0]-(1.115+0.1396))+x[0]-(1.115+0.1396)))*Polynomial(x,&par[3]);
  }

  Double_t xifit1(Double_t *x, Double_t *par) {
	  return par[0]*par[1]*Gauss(x,par[6],par[2])+
		 par[3]*par[1]*Gauss(x,par[5],par[4])+
		 Polynomial(x,&par[6]);
  }

Double_t gp1(Double_t *x, Double_t *par) {
          return 
	  par[0]*par[1]*Gauss(x,par[3],par[2])+
	  Polynomial(x,&par[3]);
}

Double_t gpExp(Double_t *x, Double_t *par) {
          return 
	  par[0]*par[1]*Gauss(x,par[3],par[2])+
		  Exp(x,&par[3]);
}



 Double_t gp2(Double_t *x, Double_t *par) {
	  return 
		  par[1]*(par[0]*par[3]/(1.+par[3])*Gauss(x,par[5],par[2])+
				 par[0]/(1.+par[3])*Gauss(x,par[5],par[4]))+
		  Polynomial(x,&par[5]);
  }

 Double_t gp22(Double_t *x, Double_t *par) {
	  return 
		  par[0]*par[1]*Gauss(x,par[3],par[2])+
		  par[0]*par[4]*Gauss(x,par[6],par[5])+
		  Exp(x,&par[6]);
  }

  Double_t g2(Double_t *x, Double_t *par) {
	  return 
		  par[1]*(par[0]*par[3]/(1.+par[3])*Gauss(x,par[5],par[2])+
			  par[0]/(1.+par[3])*Gauss(x,par[5],par[4]));
  }

  Double_t dstfit(Double_t *x, Double_t *par) {
	  return 
		  par[1]*(par[0]*par[3]/(1.+par[3])*Gauss(x,par[5],par[2])+
				 par[0]/(1.+par[3])*Gauss(x,par[5],par[4]))+
		  sqrt(0.5*(fabs(x[0]-0.1396)+x[0]-0.1396))*Polynomial(x,&par[5]);
//		  sqrt(0.5*(fabs(x[0]-0.93827231)+x[0]-0.93827231))*Polynomial(x,&par[5]);
  }


  Double_t bkg(Double_t *x, Double_t *par) {
	  double t = x[0]-0.1396;
	  double c = 2.*(0.0433333333333*par[0]+
			 0.01535976*par[1]+0.004939499188*par[2]+
			 0.00156113775127*par[3]+
			 0.000493602624474*par[4]);
	  c=1;
	  return 
		  sqrt(t)*par[0]/c*(1.+t*par[1]+
				    par[2]*pow(t,2)+
				    par[3]*pow(t,3)+
				    par[4]*pow(t,4));
		  
  }

  Double_t dstfitexp(Double_t *x, Double_t *par) {
	  return 
		  par[1]*(par[0]*par[3]/(1.+par[3])*Gauss(x,par[5],par[2])+
				 par[0]/(1.+par[3])*Gauss(x,par[6],par[4]))+
		  par[7]*(1.-exp((x[0]-0.1396)/par[8]))*pow(x[0]/0.1396,par[9]);
  }


 Double_t omegafit(Double_t *x, Double_t *par) {
	  return par[0]*par[1]*Gauss(x,par[3],par[2])+
		 sqrt(0.5*(fabs(x[0]-(1.115+0.494))+x[0]-(1.115+0.494)))*Polynomial(x,&par[3]);
  }

 Double_t xistfit(Double_t *x, Double_t *par) {
	 return par[0]*par[1]*bwg(x[0],par[3],par[2],par[4])+
		 sqrt(0.5*(fabs(x[0]-(1.321+0.1396))+x[0]-(1.321+0.1396)))*Polynomial(x,&par[4]);
 }

Double_t xist(Double_t *x, Double_t *par) {
	double dm = x[0]-(1.321+0.1396);
	if(dm<0) return 0;
	return par[0]*par[1]*bwg(x[0],par[3],par[2],par[4])+
		par[0]*pow(dm,par[5])*(par[6]+par[7]*dm+par[8]*dm*dm+
				       par[9]*dm*dm*dm);
}

Double_t xist2(Double_t *x, Double_t *par) {
	double dm = x[0]-(1.321+0.1396);
	if(dm<0) return 0;
	return par[0]*par[1]*(par[11]/(1.+par[11])*
			      bwg(x[0],par[3],par[2],par[4])+
			      bwg(x[0],par[3],par[10],par[4])/(1.+par[11]))+
		par[0]*pow(dm,par[5])*(par[6]+par[7]*dm+par[8]*dm*dm+
				       par[9]*dm*dm*dm);
}

Double_t xist1(Double_t *x, Double_t *par) {
	double dm = x[0]-(1.321+0.1396);
	if(dm<0) return 0;
	return par[0]*par[1]*(bwg(x[0],par[3],par[2],par[4])+
			      par[10]*Gauss(x,par[12],par[11]))+
		par[0]*pow(dm,par[5])*(par[6]+par[7]*dm+par[8]*dm*dm+
				       par[9]*dm*dm*dm);
}

 Double_t ddst(Double_t *x, Double_t *par) {
	 return par[0]*par[1]*bwg(x[0],par[3],par[2],par[4])+
		 par[1]*par[9]*bwg(x[0],par[10],par[2],par[11])+
		 sqrt(0.5*(fabs(x[0]-0.1396)+x[0]-0.1396))*
		  (par[5]+par[6]*(x[0]-0.1396)+par[7]*pow(x[0]-0.1396,2)+
		   par[8]*pow(x[0]-0.1396,3));
  }

 Double_t ddst1(Double_t *x, Double_t *par) {
	 Double_t dm = x[0]-(0.1396+1.869300);
	 if (dm<0) return 0;
	 return par[1]*(par[0]*bwg(x[0],par[3],par[2],par[4])+
			par[5]*Gauss(x,par[6],par[7])+
			par[8]*Gauss(x,par[9],par[10]))+
		 pow(dm,par[11])*
		 (par[12]+par[13]*dm+par[14]*dm*dm+par[15]*dm*dm*dm);
 }
//
// to fit D+ pi-
//

 Double_t dpluspi(Double_t *x, Double_t *par) {
	 Double_t dm = x[0]-(0.1396+1.869300);
	 if (dm<0) return 0;
	 return par[1]*(par[0]*bwg(x[0],par[3],par[2],par[4])+
			par[5]*Gauss(x,par[6],par[7])+
			par[8]*Gauss(x,par[9],par[10]))+
		 pow(dm,par[11])*
		 (par[12]+par[13]*dm+par[14]*dm*dm+par[15]*dm*dm*dm);
 }

//
// to fit D0 pi+
//

 Double_t d0pi(Double_t *x, Double_t *par) {
	 Double_t dm = x[0]-(0.1396+1.86450);
	 if (dm<0) return 0;
	 return par[1]*(par[0]*bwg(x[0],par[3],par[2],par[4])+
			par[5]*Gauss(x,par[6],par[7])+
			par[8]*Gauss(x,par[9],par[10]))+
		 pow(dm,par[11])*
		 (par[12]+par[13]*dm+par[14]*dm*dm+par[15]*dm*dm*dm);
 }
//
// to fit Dst pi-
//

 Double_t dstpi(Double_t *x, Double_t *par) {
	 Double_t dm = x[0]-(2.010+0.1396);
	 if (dm<0) return 0;
	 return par[1]*(par[0]*bwg(x[0],par[3],par[2],par[4])+
			par[5]*bwg(x[0],par[6],par[7],par[8]))+
		 pow(dm,par[9])*
		 (par[10]+par[11]*dm+par[12]*dm*dm+par[13]*dm*dm*dm);
 }

//
// to fit Dst p
//

 Double_t dstp(Double_t *x, Double_t *par) {
	 Double_t dm = x[0]-(2.010+0.93827231);
	 if (dm<0) return 0;
	 return par[1]*(par[0]*bwg(x[0],par[3],par[2],par[4]))+
		  sqrt(dm)*
		 (par[5]+par[6]*dm+par[7]*dm*dm+par[8]*dm*dm*dm);
 }

//
// to fit D p
//

 Double_t dp(Double_t *x, Double_t *par) {
//	 Double_t dm = x[0]-(1.86930+0.93827231);
	 Double_t dm = x[0]-(1.86450+0.93827231);
//	 Double_t dm = x[0]-(2.010+0.93827231);
	 if (dm<0) return 0;
	 return par[1]*(par[0]*bwg(x[0],par[3],par[2],par[4]))+
		 pow(dm,par[5])*
		 (par[6]+par[7]*dm+par[8]*dm*dm+par[9]*dm*dm*dm)*3.e-3;
 }

//
// to fit Dst pi-
//
 Double_t ddstbkg(Double_t *x, Double_t *par) {
	 return  sqrt(0.5*(fabs(x[0]-0.1396)+x[0]-0.1396))*
		 (par[0]+par[1]*(x[0]-0.1396)+par[2]*pow(x[0]-0.1396,2)+
		  par[3]*pow(x[0]-0.1396,3));
 }

 Double_t pentafit(Double_t *x, Double_t *par) {
	 return par[0]*par[1]*Gauss(x,par[3],par[2])+
	     sqrt(0.5*(fabs(x[0]-(1.321+0.1396))+x[0]-(1.321+0.1396)))*Polynomial(x,&par[3]);
	     
 }


 Double_t pentafit1(Double_t *x, Double_t *par) {
	 return par[0]*par[1]*bwg(x[0],par[3],par[2],par[4])+par[5]*par[1]*Gauss(x,par[7],par[6])+
		 sqrt(0.5*(fabs(x[0]-(1.321+0.1396))+x[0]-(1.321+0.1396)))*Polynomial(x,&par[7]);
	     
 }


Double_t pt(Double_t *x, Double_t *par) {
	return par[0]*exp(-x[0]/par[1])*x[0];
}

Double_t pt1(Double_t *x, Double_t *par) {
	return par[0]*exp(-x[0]/par[1])*x[0];
}

Double_t conv(Double_t *x, Double_t *par) {
	return eff(x,&par[0])*exp(par[3]*x[0]+par[4]);
}

Double_t deconv(Double_t *x, Double_t *par) {
	return exp(par[3]*x[0]+par[4])/eff(x,&par[0]);
}


Double_t lamfit(Double_t *x, Double_t *par) {
	  return 
		  par[1]*(par[0]*par[3]/(1.+par[3])*Gauss(x,par[5],par[2])+
				 par[0]/(1.+par[3])*Gauss(x,par[5],par[4]))+
		  Polynomial(x,&par[5])*sqrt(0.5*(fabs(x[0]-(0.9382723128+0.1396))+x[0]-(0.9382723128+0.1396)));
}



Double_t bethe2(Double_t *x, Double_t *par) {
	double t = x[0]/sqrt(1.+x[0]*x[0]);
	return
		(par[0]*log(fabs(x[0]/(x[0]+par[1])))+par[2])/(t*t)+par[3]*(t-1.)+par[4]*pow(t-1.,2)+par[5];
}

Double_t bethe1(Double_t *x, Double_t *par) {
	double xx = x[0]/par[5];

	double t = xx/sqrt(1.+xx*xx);
	return
		(par[0]*log(fabs(xx/(xx+par[1])))+par[2])/(t*t)+par[3]*(t-1.)+par[4]*pow(t-1.,2);
}
Double_t Gamma(Double_t *x,Double_t *par ) { 
	double m1 = par[0];
	const double m2 = 139;
	double Q=x[0];
	double M = m1+m2+Q;
	double a = m1+m2;
	double b = m1-m2;
//	double p = 10;
	double p = sqrt((M*M-a*a)*(M*M-b*b))/2./M;
//	double f = 0.75 / 93;
	double f = par[1] / 93;
	return 1./6./acos(-1.)*m1/M*f*f*p*p*p;
}

// Double_t xibfit(Double_t* x, Double_t* par){ 
// 	return 0.015*(11.136*Gauss(x,5.796,0.016)+11.910);
// }

Double_t xicfit(Double_t* x, Double_t* par){ 
	return par[0]*(par[1]*Gauss(x,par[3],par[2])+par[4]+par[5]*x[0]);
}


const double masses[] = { 
	5.799458     ,
	5.661428     ,
	5.395928     ,
	5.793429     ,
	5.787025     ,
	5.813913     ,
	5.801238     ,
	5.796516     ,
	5.795753     ,
	5.790754     ,
	5.807012     ,
	5.923625     ,
	5.778005     ,
	5.756891     ,
	6.325589     ,
	5.800344     ,
	5.872428     ,
	5.774071     ,
	6.223488     ,
	5.791504     ,
	5.765860     ,
	5.816276     ,
	5.886657     ,
	5.796786     ,
	5.253847     ,
	5.781273 };

const double errors[] = {
	0.2508454E-02 ,
	0.9936395E-02 ,
	0.5969557E-02 ,
	0.3743158E-02 ,
	0.5364847E-02 ,
	0.9280185E-02 ,
	0.5981869E-02 ,
	0.8329995E-02 ,
	0.7055258E-02 ,
	0.9665354E-02 ,
	0.1182285E-01 ,
	0.4683272E-02 ,
	0.6402510E-02 ,
	0.1055388E-01 ,
	0.5321487E-02 ,
	0.8818955E-02 ,
	0.9964406E-02 ,
	0.7545226E-02 ,
	0.5445669E-02 ,
	0.2266543E-02 ,
	0.5950275E-02 ,
	0.8225833E-02 ,
	0.1155143E-01,
	0.4386774E-02 ,
	0.5132507E-02 ,
	0.6774143E-02 };


Double_t xibfit(Double_t* x, Double_t* par){ 
	double sum=0.;
	for (int i=0;i<26;i++) { 
		sum+=Gauss(x,par[3],par[2]*errors[i]);
	}
	sum/=26.;
	return par[0]*(par[1]*sum+par[4]+par[5]*x[0]);
}


Double_t xib2gfit(Double_t* x, Double_t* par){ 
	return par[0]*(par[1]*(Gauss(x,par[4],par[2])*par[5]/(1.+par[5])+
			       Gauss(x,par[4],par[3])/(1.+par[5]))+
		       par[6]+par[7]*x[0]);
}


Double_t line(Double_t* x, Double_t* par){ 
	static const double Q=5796.0-3096.916-1321.34;
	return par[1]*(x[0]-Q)+par[0];
}

