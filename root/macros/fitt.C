{ 
        gInterpreter->ProcessLine(".L fit.C");
        gStyle->SetOptFit(111);
	TH1F* h = new TH1F("dstpi","dstpi",70,2.21,2.91);
	const int bins[] = { 222, 266, 276, 259, 265, 277, 289, 269, 294, 264, 276, 288, 280, 289, 285, 276, 274, 286, 337, 349, 348, 390, 321, 332, 312, 331, 297, 282, 275, 272, 273, 263, 272, 259, 260, 251, 261, 245, 233, 240, 239, 226, 228, 247, 234, 218, 212, 218, 219, 215, 207, 197, 213, 211, 220, 204, 198, 204, 189, 202, 195, 186, 215, 191, 197, 187, 172, 185, 190, 171 };
	int nbins = sizeof(bins)/sizeof(int);
	for (int i=1; i<=h->GetNbinsX(); i++) { 
		h->SetBinContent(i,bins[i-1]);
	}
	
	TF1* fitFun = new TF1("gp",FitFunction, h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),11);

	fitFun->SetParameter(0,0);
	fitFun->SetParameter(1,10);
	fitFun->SetParameter(2,0);
	fitFun->SetParameter(3,2420);
	fitFun->SetParameter(4,0);
	fitFun->SetParameter(5,10);
	fitFun->SetParameter(6,0);
	fitFun->SetParameter(7,2460);
	fitFun->SetParLimits(0,1,1);
	fitFun->SetParLimits(1,1,1);
	fitFun->SetParLimits(2,1,1);
	fitFun->SetParLimits(3,1,1);
	fitFun->SetParLimits(4,1,1);
	fitFun->SetParLimits(5,1,1);
	fitFun->SetParLimits(6,1,1);
	fitFun->SetParLimits(7,1,1);
	h->Fit("gp","M");
        fitFun->SetParameter(2,18.9);
        fitFun->SetParameter(6,23);
	h->Fit("gp","MV");
	fitFun->SetParLimits(0,0,1000);
	fitFun->SetParLimits(4,0,1000);
	h->Fit("gp","MV");
	fitFun->SetParLimits(3,2410,2440);
	fitFun->SetParLimits(7,2445,2470);
	h->Fit("gp","MV");
	

        gPad->SaveAs("root.ps");


}

