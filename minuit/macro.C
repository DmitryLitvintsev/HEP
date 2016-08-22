{
   TH1F* h  = new TH1F("width","width",100,-10,10);
   TH1F* h1 = new TH1F("mass","mass",100,-10,10);

   FILE* fp = fopen("residuals_fixed_sigma_10_gausian_sigma_range.dat","r");
   float a,b,c;
   {while (fscanf(fp,"%f %f %f",&a,&b,&c)!=EOF) { 
	   h->Fill(a);
	   h1->Fill(b);
   }}

   TH1F h100(*h);
   TH1F h101(*h1);
   h100.SetName("width1");
   h101.SetName("mass1");

   TH1F* h  = new TH1F("width","width",100,-10,10);
   TH1F* h1 = new TH1F("mass","mass",100,-10,10);

   FILE* fp = fopen("residuals_fixed_width_3.2_1.4-2.2.dat","r");
   float a,b,c;
   {while (fscanf(fp,"%f %f %f",&a,&b,&c)!=EOF) { 
	   h->Fill(a);
	   h1->Fill(b);
   }}

   TH1F h200(*h);
   TH1F h201(*h1);
   h200.SetName("width2");
   h201.SetName("mass2");


   TH1F* h  = new TH1F("width","width",100,-10,10);
   TH1F* h1 = new TH1F("mass","mass",100,-10,10);

   FILE* fp = fopen("residuals_fixed_sigma_0-3.dat","r");
   float a,b,c;
   {while (fscanf(fp,"%f %f %f",&a,&b,&c)!=EOF) { 
	   h->Fill(a);
	   h1->Fill(b);
   }}

   TH1F h300(*h);
   TH1F h301(*h1);
   h300.SetName("width3");
   h301.SetName("mass3");


   TH1F* h  = new TH1F("a","a",100,-10,10);
   TH1F* h1 = new TH1F("w","w",100,-10,10);
   TH1F* h2 = new TH1F("g","g",100,-10,10);
   TH1F* h3 = new TH1F("m","m",100,-10,10);
   TH1F* h4 = new TH1F("pwr","pwr",100,-10,10);

   TH1F* h5 = new TH1F("a0","a0",100,-100,100);
   TH1F* h6 = new TH1F("a1","a1",100,-100,100);
   TH1F* h7 = new TH1F("a2","a2",100,-100,100);
   TH1F* h8 = new TH1F("a3","a3",100,-100,100);

   FILE* fp = fopen("toy_residuals.txt","r");
   float a,w,g,m,pwr,a0,a1,a2,a3;
   {while (fscanf(fp,"%f %f %f %f %f %f %f %f %f",&a,&w,&g,&m,&pwr,&a0,&a1,&a2,&a3)!=EOF) { 
	   h->Fill(a);
	   h1->Fill(w);
	   h2->Fill(g);
	   h3->Fill(m);
	   h4->Fill(pwr);
	   h5->Fill(a0);
	   h6->Fill(a1);
	   h7->Fill(a2);
	   h8->Fill(a3);
   }}

   

}
x
