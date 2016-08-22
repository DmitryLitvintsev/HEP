void draw_polyline(TH1*h) { 
	TPolyLine* line = new  TPolyLine(h->GetNbinsX());
	for (int i=0; i<h->GetNbinsX(); i++) { 
		double y = h->GetBinContent(i+1);
		double x = h->GetBinCenter(i+1);
		line->SetPoint(i,x,y);
	}
	line->SetLineColor(2);
	line->SetLineWidth(2);
	line->Draw("same");
}
