Double_t bethe(Double_t *x, Double_t *par) {
	double t = x[0]/sqrt(1.+x[0]*x[0]);
	return
		(par[0]*log(fabs(x[0]/(x[0]+par[1])))+par[2])/(t*t)+par[3]*(t-1.)+par[4]*pow(t-1.,2);
}
