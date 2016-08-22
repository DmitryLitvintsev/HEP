
int add(TH1F* one, TH1F* two, TH1F* result, double a, double b) { 
	if (one->GetNbinsX()!=two->GetNbinsX()) return NULL;
	for (int i=1;i<=one->GetNbinsX(); i++) { 
		double value = a*one->GetBinContent(i)+b*two->GetBinContent(i);
		double error = sqrt(a*a*one->GetBinContent(i)+b*b*two->GetBinContent(i));
		result->SetBinContent(i,value);
		result->SetBinError(i,error);
		
	}
	return 0;
}
