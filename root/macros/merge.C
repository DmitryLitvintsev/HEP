{
	TChain* xi  = new TChain("Xic/hyperon");
        xi->Add("h*.root");
	TTree* newtree = xi->CopyTree("","");
	TFile f("sum_1.root","NEW");
	f.cd();
	newtree->Write();
	f.Write();
	f.Close();
}
