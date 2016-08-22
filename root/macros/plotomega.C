int plotomega(TH1*h) {
	TFrame* frame = (TFrame*) gPad->GetListOfPrimitives()->FindObject("TFrame");
	double x1 = 1.775;
	double x2 = 1.8;
	double y = (frame->GetY2()-frame->GetY1());
	double y12 = frame->GetY2()-0.1*y;
	double height = 0.08*y;
	double y11 = y12-height;
	double delta = 0.3*height;
	double y22 = y11 - delta;
	double y21 = y22 - height;

	TPave* p1 = new TPave(x1,y11,x2,y12);
	TPave* p2 = new TPave(x1,y21,x2,y22);
	
	p1->SetBorderSize(1);
	p1->SetFillColor(0);
	p2->SetBorderSize(1);
	p2->SetFillColor(5);
	p2->SetLineColor(5);
	
	TLatex* txt1 = new TLatex(x2+0.0075,y12-0.5*height,"#Lambda K^{-}");
	TLatex* txt2 = new TLatex(x2+0.0075,y22-0.5*height,"#Lambda K^{+}");
	
	txt1->SetTextAlign(12);
	txt2->SetTextAlign(12);

	p1->Draw();
	p2->Draw();
	txt1->Draw();
	txt2->Draw();
	return 0;
	
}
