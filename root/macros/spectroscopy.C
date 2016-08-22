{
  c1 = new TCanvas("c1","spectroscopy",10,10,500,600);
  c1->Range(-1,0,1.2,0.97);

  TGaxis *axis = new TGaxis(-0.5,0.05,-0.5,0.95,2.2,3.0,20210,"");

  double offset = 0.05;
  double scale = (3.0 - 2.2 ) / (0.95 - 0.05);

  axis->SetLabelOffset(0.01);
  axis->SetLabelSize(0.05);
  axis->SetName("axis");
  axis->SetTitle("mass [GeV]             ");
  axis->SetTitleSize(0.08);
  axis->SetTitleOffset(1.2);
  axis->SetTextFont(42);
  axis->SetLabelFont(42);
  axis->Draw();

  TLatex* t1 = new TLatex(-0.3,0.01,"#Lambda_{c} - type"); 
  t1->SetTextFont(42);
  t1->SetTextAlign(11);
  t1->Draw();

  TLatex* t2 = new TLatex(0.4,0.01,"#Sigma_{c} - type"); 
  t2->SetTextFont(42);
  t2->SetTextAlign(11);
  t2->Draw();

  double lw = 0.4;
  double xstart = -0.3;

  double y   = offset + ( 2.2849 - 2.2) / scale;
  double x1  = xstart; 
  double x2  = xstart + lw;

  TLine* l0 = new TLine(x1,y,x2,y);
  l0->SetLineWidth(2);
  l0->Draw();

   TLatex* txt0 = new TLatex(x2,y,"1/2^{+}");
   txt0->SetTextAlign(12);
   txt0->SetTextFont(42);
   txt0->Draw();
   txt0->SetTextSize(0.03);

  y   = offset + ( 2.5939 - 2.2) / scale;
  x2  = xstart + lw;

  TLine* l1 = new TLine(x1,y,x2,y);
  l1->SetLineWidth(2);
  l1->Draw();

   TLatex* txt0 = new TLatex(x2,y,"1/2^{-}");
   txt0->SetTextAlign(12);
   txt0->SetTextFont(42);
   txt0->Draw();
   txt0->SetTextSize(0.03);


  double aoffset = lw/8.;
  double asize = 0.01;

   TArrow* a1 = new TArrow(l1->GetX1()+aoffset,l1->GetY1(),l0->GetX1()+aoffset,l0->GetY2(),asize);
   a1->Draw();

   TArrow* a2 = new TArrow(l1->GetX1()+2.*aoffset,l1->GetY1(),l0->GetX1()+2.*aoffset,l0->GetY2(),asize);
   a2->Draw();

  y   = offset + ( 2.6266 - 2.2) / scale;
  x2  = xstart + lw;

   TLine* l2 = new TLine(x1,y,x2,y);
   l2->SetLineWidth(2);
   l2->Draw();

   TLatex* txt0 = new TLatex(x2,y,"3/2^{-}");
   txt0->SetTextAlign(12);
   txt0->SetTextFont(42);
   txt0->Draw();
   txt0->SetTextSize(0.03);



   TArrow* a3 = new TArrow(l2->GetX1()+4.*aoffset,l2->GetY1(),l0->GetX1()+4.*aoffset,l0->GetY2(),asize);
   a3->Draw();

   TArrow* a4 = new TArrow(l2->GetX1()+5.*aoffset,l2->GetY1(),l0->GetX1()+5.*aoffset,l0->GetY2(),asize);
   a4->Draw();

   TLatex* twopi = new TLatex(-0.15,offset+(0.5*(2.626+2.285)-2.2)/scale,"2#pi");
   twopi->SetTextFont(42);
   twopi->SetTextAlign(22);
   twopi->Draw();

   TLatex* twopi = new TLatex(-0.15,offset+(0.5*(2.626+2.285)-2.2)/scale,"2#pi");
   twopi->SetTextFont(42);
   twopi->SetTextAlign(22);
   twopi->Draw();

   TLatex* dwave = new TLatex(0.28,0.45,"D-wave");
   dwave->SetTextFont(42);
   dwave->SetTextAlign(11);
   dwave->Draw();

   TLatex* swave = new TLatex(0.0,0.37,"S-wave");
   swave->SetTextFont(42);
   swave->SetTextAlign(11);
   swave->Draw();

   TLatex* pwave = new TLatex(0.32,0.22,"P-wave");
   pwave->SetTextFont(42);
   pwave->SetTextAlign(11);
   pwave->Draw();



//----------------------------------------------------------------------

  xstart = 0.4;
  x1 = xstart;


  y   = offset + ( 2.452 - 2.2) / scale;
  x2  = xstart + lw;

  TLine* l3 = new TLine(x1,y,x2,y);
  l3->SetLineWidth(2);
  l3->Draw();

   TLatex* txt0 = new TLatex(x2,y,"1/2^{+}");
   txt0->SetTextAlign(12);
   txt0->SetTextFont(42);
   txt0->Draw();
   txt0->SetTextSize(0.03);

  y   = offset + ( 2.518 - 2.2) / scale;
  x2  = xstart + lw;

  TLine* l4 = new TLine(x1,y,x2,y);
  l4->SetLineWidth(2);
  l4->Draw();

   TLatex* txt0 = new TLatex(x2,y,"3/2^{+}");
   txt0->SetTextAlign(12);
   txt0->SetTextFont(42);
   txt0->Draw();
   txt0->SetTextSize(0.03);

   TArrow* a5 = new TArrow(l2->GetX1()+7.*aoffset,l2->GetY1(),l3->GetX1()+2.8*aoffset,l3->GetY2(),asize);
   a5->Draw();

   TArrow* a4 = new TArrow(l1->GetX1()+7.*aoffset,l1->GetY1(),l3->GetX1()+aoffset,l3->GetY2(),asize);
   a4->Draw();


   TArrow* a6 = new TArrow(l3->GetX1()+aoffset,l3->GetY1(),l0->GetX1()+6.*aoffset,l0->GetY2(),asize);
   a6->Draw();

   TArrow* a6 = new TArrow(l4->GetX1()+5.2*aoffset,l4->GetY1(),l0->GetX1()+7.*aoffset,l0->GetY2(),asize);
   a6->Draw();

   y   = offset + ( 2.765 - 2.2) / scale;
   x2  = xstart + lw;

   TLatex* txt0 = new TLatex(x2,y-0.002,"1/2^{-}");
   txt0->SetTextAlign(13);
   txt0->SetTextFont(42);
   txt0->Draw();
   txt0->SetTextSize(0.03);

   TLine* l5 = new TLine(x1,y,x2,y);
   l5->SetLineWidth(2);
   l5->Draw();

   TLatex* txt0 = new TLatex(x2,y+0.002,"1/2^{-},3/2^{-}");
   txt0->SetTextAlign(11);
   txt0->SetTextFont(42);
   txt0->Draw();
   txt0->SetTextSize(0.03);

   y   = offset + ( 2.770 - 2.2) / scale;
   x2  = xstart + lw;

   TLine* l6 = new TLine(x1,y,x2,y);
   l6->SetLineWidth(2);
   l6->Draw();


   y   = offset + ( 2.805 - 2.2) / scale;
   x2  = xstart + lw;

   TLine* l7 = new TLine(x1,y,x2,y);
   l7->SetLineWidth(2);
   l7->Draw();

   TLatex* txt0 = new TLatex(x2,y,"3/2^{-}");
   txt0->SetTextAlign(13);
   txt0->SetTextFont(42);
   txt0->Draw();
   txt0->SetTextSize(0.03);

   y   = offset + ( 2.815 - 2.2) / scale;
   x2  = xstart + lw;

   TLine* l8 = new TLine(x1,y,x2,y);
   l8->SetLineWidth(2);
   l8->Draw();

   TLatex* txt0 = new TLatex(x2,y,"5/2^{-}");
   txt0->SetTextAlign(11);
   txt0->SetTextFont(42);
   txt0->Draw();
   txt0->SetTextSize(0.03);


  

}
