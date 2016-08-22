TF1* f = new TF1("Gamma",Gamma,0,100,2);
f->SetParameter(0,5619);
f->SetParameter(1,0.75);
f->Draw();

TF1* f1 = new TF1("Gamma1",Gamma,0,100,1);
f1->SetParameter(0,2285);
f1->Draw();

f->Draw();
f1->Draw("same");

float m_pi = 139.57018;

float x[]   = { 167.587-139.57018, 234.5-139.57018};
float y[]   = { 2.23, 17.9 };
float dy[]  = { 0.3 , sqrt(3.8*3.8+4.*4.) };

TGraphErrors* gr = new TGraphErrors(2,x,y,0,dy);

//gr->SetMarkerColor(4);
gr->SetMarkerStyle(20);
gr->Draw("AP");


