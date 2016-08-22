int plotxi1(TH1*h) {
	TFrame* frame = (TFrame*) gPad->GetListOfPrimitives()->FindObject("TFrame");
//	if (gPad->GetLogy()) { 
// 	double x1 = 1.34;
// 	double x2 = 1.36;
// 	double x1 = 2.6;
// 	double x2 = 2.65;
// 	double x1 = 1.82;
// 	double x2 = 1.84;
// 	double x1 = 4.5;
// 	double x2 = 4.9;
//   	double x1 = 1.8;
//   	double x2 = 1.85;

// 		double x1 = 10;
// 		double x2 = 12;
// 		double y = (exp(frame->GetY2())-exp(frame->GetY1()));
// 		double y12 = exp(frame->GetY2())-0.05*y;
// 		double height = 0.08*y;
// 		double y11 = y12-height;
// 		double delta = 0.3*height;
// 		double y22 = y11 - delta;
// 		double y21 = y22 - height;
		
// 		TPave* p1 = new TPave(x1,y11,x2,y12);
// 		TPave* p2 = new TPave(x1,y21,x2,y22);

// 		p1->SetBorderSize(1);
// 		p1->SetFillColor(2);
// 		p1->SetLineWidth(2);
// 		p2->SetLineWidth(2);
// 		p2->SetBorderSize(1);
// 		p2->SetFillColor(5);
// 		p2->SetLineColor(1);
		
// 		TLatex* txt1 = new TLatex(x2+0.1,y22-0.5*height,"#Lambda#pi^{+}");
// 		TLatex* txt2 = new TLatex(x2+0.1,y12-0.5*height,"#Lambda#pi^{-}");
// 		txt1->SetTextAlign(12);
// 		txt2->SetTextAlign(12);

// 		p1->Draw();
// 		p2->Draw();
// 		txt1->Draw();
// 		txt2->Draw();
// 	}
// 	else { 
		
		double x1 = 10;
		double x2 = 12;
		double y = (frame->GetY2()-frame->GetY1());
		double y12 = frame->GetY2()-0.05*y;
		double height = 0.08*y;
		double y11 = y12-height;
		double delta = 0.3*height;
		double y22 = y11 - delta;
		double y21 = y22 - height;

		TPave* p1 = new TPave(x1,y11,x2,y12);
		TPave* p2 = new TPave(x1,y21,x2,y22);



		
		p1->SetBorderSize(1);
		p1->SetFillColor(1);
		p1->SetFillStyle(3001);
		p1->SetLineWidth(2);
		p2->SetLineWidth(2);
		p2->SetBorderSize(1);
		p2->SetFillColor(5);
		p2->SetLineColor(1);
	

		TLatex* txt1 = new TLatex(x2+0.1,y22-0.5*height,"#Lambda#pi^{+}");
		TLatex* txt2 = new TLatex(x2+0.1,y12-0.5*height,"#Lambda#pi^{-}");

		cout << " hi " << endl;
		
		txt1->SetNDC();
		txt2->SetNDC();
		
		txt1->SetTextAlign(12);
		txt2->SetTextAlign(12);
		
		p1->Draw();
		p2->Draw();
		txt1->Draw();
		txt2->Draw();
//	}
	return 0;
	
}
