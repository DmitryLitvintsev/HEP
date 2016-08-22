//
// COT dEdx offline corrections & UC (root version)
//
// Version 4.0 - July 13th, 2004:
//   - CDF6932 track phi0, eta, usedhits corr.
//   - CDFXXXX track momentum dependent corrections
//   - Fix for the COT pressure corr. for runs>160457 (5.X data)
//   - Universal curves: XXXX
//   - dEdx Resolution model (includes mometum dependence): XXXX
//
// D.Ambrose, P.Catastini, M.Donega, S.Giagu, M.Jones, J.D.Lewis, 
// N.Pounder, G.Punzi, P.Squillacioti, V.Tiwari, D.Tonelli, S.S.Eiko Yu
//

//Get Mass given PDG ID
Double_t getMass(Int_t pid) {
  Double_t mass=-999;
  switch(abs(pid)) { 
    case 2212: mass= 0.938272; break;// proton
    case 211:  mass= 0.139570; break;// pion
    case 321:  mass= 0.493677; break;// kaon
    case 13:   mass= 0.105658; break;// muon
    case 11:   mass= 0.000511; break;// electron
    default: mass = -999;	
  }	
  return mass;
}

//Return Run bin index (used for time dependent corrections)
Int_t RunIndex(Int_t run) {

  Int_t runRange[11] = {142558, 146200, 150000, 
		      152635, 154900, 160458, 
		      161950, 164275, 165900, 
		      167830, 168889};

  for (Int_t j=0; j<11; j++) {
    if (run <= runRange[j]) {
      return j;
    }
  }
  return 10;
}

//Retrun Phi bin index (used for phi0 corrections)
Int_t PhiIndex(Double_t phi) {
  Int_t indexphi = static_cast<Int_t>(phi*180./3.141592654/4.);
  if (indexphi < 0)  indexphi = 0;
  if (indexphi > 89) indexphi = 89;
  return indexphi;
}

//Eta Correction (postive charge)
//V3.0 02/24/2004 SG 
Double_t EtaCorrPos(Int_t run, Double_t eta) {

  //Time independent correction
  // corr = 15/(eta_c[0] + eta_c[1]*eta + ... + eta_c[8]*eta^8)

  Double_t eta_c[9] = {13.010, -0.66345, 4.8624, 1.4566, -19.520, 
		     -2.9009, 30.600, 1.8635, -15.204};
  
  Double_t corr1  = 1.;
  Double_t eta_tr = eta;

  if (eta_tr >  1.) eta_tr =  1.; //Diego's corrections def. in [-1,1]
  if (eta_tr < -1.) eta_tr = -1.;

  Double_t invcorr1 = eta_c[0] + eta_c[1]*eta_tr;
  for (Int_t j=2; j<9; j++) 
    invcorr1 += eta_c[j]*pow(eta_tr, j);
  if (invcorr1 != 0) corr1 = 15./invcorr1;
  
  //Time dependent correction
  // Corr = 15/(f1(run index) + f2(run index)*eta + f3(run index)*eta^2)
  
  Double_t escale[3][11]; //3 pol. coeff. * 11 run ranges

  escale[0][0]  =  14.989;
  escale[0][1]  =  14.516;
  escale[0][2]  =  15.292;
  escale[0][3]  =  15.026;
  escale[0][4]  =  15.014;
  escale[0][5]  =  15.020;
  escale[0][6]  =  14.959;
  escale[0][7]  =  14.757;
  escale[0][8]  =  14.950;
  escale[0][9]  =  14.922;
  escale[0][10] =  15.116;
                             
  escale[1][0]  =  0.16654;
  escale[1][1]  =  0.10031;
  escale[1][2]  =  0.68078;
  escale[1][3]  =  0.45194;
  escale[1][4]  =  0.58307;
  escale[1][5]  =  0.17333;
  escale[1][6]  = -0.48656e-01;
  escale[1][7]  = -0.10283;
  escale[1][8]  = -0.28809;
  escale[1][9]  = -0.45242;
  escale[1][10] = -0.25658;

  escale[2][0]  =  0.;
  escale[2][1]  =  2.4964;
  escale[2][2]  =  0.17934;
  escale[2][3]  =  0.29188;
  escale[2][4]  =  0.14810;
  escale[2][5]  =  0.;
  escale[2][6]  =  0.;
  escale[2][7]  =  0.;
  escale[2][8]  =  0.;
  escale[2][9]  =  0.;
  escale[2][10] = -0.73837;

  Int_t runIndex = RunIndex(run);

  Double_t corr2    = 1.;
  Double_t invcorr2 = escale[0][runIndex] + escale[1][runIndex]*eta_tr + 
                    escale[2][runIndex]*eta_tr*eta_tr;
  if (invcorr2 != 0) corr2 = 15./invcorr2;
    
  
  //Total eta correction
  Double_t totcorr = corr1*corr2;

  return totcorr;
}

//Eta Correction (negative charge)
//V3.0 02/24/2004 SG 
Double_t EtaCorrNeg(Int_t run, Double_t eta) {

  //Time independent correction
  // corr = 15/(eta_c[0] + eta_c[1]*eta + ... + eta_c[9]*eta^9)

  Double_t eta_c[10] = {12.958, -0.77330, 4.5574, 1.5125, -17.331, 
                      3.0136, 27.124, -13.513, -13.619, 9.2116};

  Double_t corr1    = 1.;
  Double_t eta_tr = eta;

  if (eta_tr >  1.) eta_tr =  1.; //Diego's corrections def. in [-1,1]
  if (eta_tr < -1.) eta_tr = -1.;

  Double_t invcorr1 = eta_c[0] + eta_c[1]*eta_tr;
  for (Int_t j=2; j<10; j++) 
    invcorr1 += eta_c[j]*pow(eta_tr, j);
  if (invcorr1 != 0) corr1 = 15./invcorr1;
  
  
  //Time dependent correction
  // Corr = 15/(f1(run index) + f2(run index)*eta + f3(run index)*eta^2)
  
  Double_t escale[3][11]; //3 pol. coeff. * 11 run ranges

  escale[0][0]  = 15.160;
  escale[0][1]  = 14.885;
  escale[0][2]  = 15.121;
  escale[0][3]  = 15.143;
  escale[0][4]  = 14.961;
  escale[0][5]  = 14.904;
  escale[0][6]  = 14.896;
  escale[0][7]  = 14.944;
  escale[0][8]  = 14.966;
  escale[0][9]  = 15.000;
  escale[0][10] = 15.110;
  
  escale[1][0]  = 0.64821;
  escale[1][1]  = -0.56363e-03;
  escale[1][2]  = 0.36496;
  escale[1][3]  = 0.61759;
  escale[1][4]  = 0.50315;
  escale[1][5]  = 0.30209;
  escale[1][6]  = 0.038832;
  escale[1][7]  = -0.10081;
  escale[1][8]  = -0.20159;
  escale[1][9]  = -0.54208;
  escale[1][10] = -0.55672 ;
  
  escale[2][0]  = 0.;
  escale[2][1]  = 0.;
  escale[2][2]  = 0.;
  escale[2][3]  = 0.;
  escale[2][4]  = 0.;
  escale[2][5]  = 0.;
  escale[2][6]  = 0.;
  escale[2][7]  = 0.;
  escale[2][8]  = 0.;
  escale[2][9]  = -0.40844;
  escale[2][10] = -1.3372;

  Int_t runIndex = RunIndex(run);

  Double_t corr2    = 1.;
  Double_t invcorr2 = escale[0][runIndex] + escale[1][runIndex]*eta_tr + 
                    escale[2][runIndex]*eta_tr*eta_tr;
  if (invcorr2 != 0) corr2 = 15./invcorr2;
    
  
  //Total eta correction
  Double_t totcorr = corr1*corr2;

  return totcorr;
}

//Correction for the number of used hits (positive charge)
//V3.0 02/24/2004 SG 
Double_t HitsCorrPos(Int_t nhit) {

  // corr = 15/(hit_c[0] + hit_c[1]*nhit + ... + eta_c[4]*nhit^4)

  Double_t hit_c[5] = {95.765, -5.0503, 0.12254, -0.13379e-02, 0.547e-05};

  Int_t nhit_tr = nhit;
  if (nhit_tr < 43) nhit_tr = 43; //Diego's corr. def. in [43,...]

  Double_t corr    = 1.;
  Double_t invcorr = hit_c[0] + hit_c[1]*static_cast<Double_t>(nhit_tr);
  for (Int_t j=2; j<5; j++) 
    invcorr += hit_c[j]*pow(static_cast<Double_t>(nhit_tr), j);
  if (invcorr != 0) corr = 15./invcorr;

  return corr;
}

//Correction for the number of used hits (negative charge)
//V3.0 02/24/2004 SG 
Double_t HitsCorrNeg(Int_t nhit) {

  // corr = 15/(hit_c[0] + hit_c[1]*nhit + eta_c[2]*nhit^2)

  Double_t hit_c[3] = {21.994, -0.14710, 0.00069498};

  Int_t nhit_tr = nhit;
  if (nhit_tr < 43) nhit_tr = 43; //Diego's corr. def. in [43,...]

  Double_t corr    = 1.;
  Double_t invcorr = hit_c[0] + hit_c[1]*nhit_tr + hit_c[2]*nhit_tr*nhit_tr;
  if (invcorr != 0) corr = 15./invcorr;

  return corr;
}

//Phi Correction (positive charge)
//V3.0 02/24/2004 SG 
Double_t PhiCorrPos(Int_t run, Double_t phi) {

  // corr = phiscale(index phi, index run)

  Double_t phiscale[90][11]; //90 bins in phi * 11 run ranges
  phiscale[0][0] = 1.02246;
  phiscale[1][0] = 0.987088;
  phiscale[2][0] = 0.992427;
  phiscale[3][0] = 0.982566;
  phiscale[4][0] = 1.0034;
  phiscale[5][0] = 0.966922;
  phiscale[6][0] = 0.987228;
  phiscale[7][0] = 0.996939;
  phiscale[8][0] = 0.953668;
  phiscale[9][0] = 0.967775;
  phiscale[10][0] = 0.958193;
  phiscale[11][0] = 0.951547;
  phiscale[12][0] = 0.988446;
  phiscale[13][0] = 0.938814;
  phiscale[14][0] = 0.951362;
  phiscale[15][0] = 0.962854;
  phiscale[16][0] = 0.945072;
  phiscale[17][0] = 0.972096;
  phiscale[18][0] = 1.00073;
  phiscale[19][0] = 0.933906;
  phiscale[20][0] = 0.96581;
  phiscale[21][0] = 0.951805;
  phiscale[22][0] = 0.940652;
  phiscale[23][0] = 0.93332;
  phiscale[24][0] = 0.926;
  phiscale[25][0] = 0.913541;
  phiscale[26][0] = 0.910077;
  phiscale[27][0] = 0.897615;
  phiscale[28][0] = 0.915909;
  phiscale[29][0] = 0.92432;
  phiscale[30][0] = 0.94555;
  phiscale[31][0] = 0.937281;
  phiscale[32][0] = 0.937643;
  phiscale[33][0] = 0.924013;
  phiscale[34][0] = 0.935159;
  phiscale[35][0] = 0.901979;
  phiscale[36][0] = 0.93761;
  phiscale[37][0] = 0.972106;
  phiscale[38][0] = 0.988773;
  phiscale[39][0] = 0.967586;
  phiscale[40][0] = 0.984265;
  phiscale[41][0] = 0.93835;
  phiscale[42][0] = 0.945438;
  phiscale[43][0] = 0.939731;
  phiscale[44][0] = 0.968187;
  phiscale[45][0] = 0.981587;
  phiscale[46][0] = 0.95407;
  phiscale[47][0] = 0.972408;
  phiscale[48][0] = 0.965809;
  phiscale[49][0] = 0.966917;
  phiscale[50][0] = 0.992609;
  phiscale[51][0] = 0.975534;
  phiscale[52][0] = 0.959696;
  phiscale[53][0] = 0.947515;
  phiscale[54][0] = 0.978949;
  phiscale[55][0] = 0.961308;
  phiscale[56][0] = 0.924827;
  phiscale[57][0] = 0.955296;
  phiscale[58][0] = 0.976715;
  phiscale[59][0] = 0.975148;
  phiscale[60][0] = 0.966676;
  phiscale[61][0] = 0.973533;
  phiscale[62][0] = 0.981655;
  phiscale[63][0] = 0.94459;
  phiscale[64][0] = 0.951827;
  phiscale[65][0] = 0.974544;
  phiscale[66][0] = 0.968989;
  phiscale[67][0] = 1.00552;
  phiscale[68][0] = 0.980609;
  phiscale[69][0] = 0.968208;
  phiscale[70][0] = 0.972445;
  phiscale[71][0] = 0.990251;
  phiscale[72][0] = 1.02835;
  phiscale[73][0] = 1.03123;
  phiscale[74][0] = 0.996117;
  phiscale[75][0] = 0.984166;
  phiscale[76][0] = 0.980236;
  phiscale[77][0] = 0.991871;
  phiscale[78][0] = 0.977526;
  phiscale[79][0] = 0.975209;
  phiscale[80][0] = 0.955258;
  phiscale[81][0] = 0.971304;
  phiscale[82][0] = 1.01284;
  phiscale[83][0] = 0.967585;
  phiscale[84][0] = 0.980261;
  phiscale[85][0] = 0.941046;
  phiscale[86][0] = 0.973774;
  phiscale[87][0] = 0.989351;
  phiscale[88][0] = 0.957441;
  phiscale[89][0] = 1.00544;
  phiscale[0][1] = 0.963792;
  phiscale[1][1] = 0.957252;
  phiscale[2][1] = 0.924546;
  phiscale[3][1] = 0.951566;
  phiscale[4][1] = 0.959669;
  phiscale[5][1] = 0.977291;
  phiscale[6][1] = 0.939379;
  phiscale[7][1] = 0.933459;
  phiscale[8][1] = 0.933524;
  phiscale[9][1] = 0.917699;
  phiscale[10][1] = 0.923302;
  phiscale[11][1] = 0.950908;
  phiscale[12][1] = 0.935824;
  phiscale[13][1] = 0.923316;
  phiscale[14][1] = 0.936158;
  phiscale[15][1] = 0.942281;
  phiscale[16][1] = 0.931186;
  phiscale[17][1] = 0.934478;
  phiscale[18][1] = 0.941474;
  phiscale[19][1] = 0.90997;
  phiscale[20][1] = 0.959283;
  phiscale[21][1] = 0.948634;
  phiscale[22][1] = 0.907511;
  phiscale[23][1] = 0.914463;
  phiscale[24][1] = 0.876387;
  phiscale[25][1] = 0.892637;
  phiscale[26][1] = 0.899999;
  phiscale[27][1] = 0.878808;
  phiscale[28][1] = 0.901514;
  phiscale[29][1] = 0.910708;
  phiscale[30][1] = 0.92175;
  phiscale[31][1] = 0.887335;
  phiscale[32][1] = 0.906839;
  phiscale[33][1] = 0.89514;
  phiscale[34][1] = 0.926312;
  phiscale[35][1] = 0.888332;
  phiscale[36][1] = 0.931757;
  phiscale[37][1] = 0.928202;
  phiscale[38][1] = 0.907693;
  phiscale[39][1] = 0.945831;
  phiscale[40][1] = 0.94983;
  phiscale[41][1] = 0.944133;
  phiscale[42][1] = 0.900333;
  phiscale[43][1] = 0.909877;
  phiscale[44][1] = 0.902951;
  phiscale[45][1] = 0.94595;
  phiscale[46][1] = 0.918038;
  phiscale[47][1] = 0.922418;
  phiscale[48][1] = 0.89282;
  phiscale[49][1] = 0.971695;
  phiscale[50][1] = 0.980133;
  phiscale[51][1] = 0.962145;
  phiscale[52][1] = 0.951315;
  phiscale[53][1] = 0.908982;
  phiscale[54][1] = 0.941874;
  phiscale[55][1] = 0.945233;
  phiscale[56][1] = 1.00449;
  phiscale[57][1] = 0.952351;
  phiscale[58][1] = 0.926517;
  phiscale[59][1] = 0.978972;
  phiscale[60][1] = 0.90433;
  phiscale[61][1] = 0.936911;
  phiscale[62][1] = 0.952196;
  phiscale[63][1] = 0.916912;
  phiscale[64][1] = 0.923133;
  phiscale[65][1] = 0.953089;
  phiscale[66][1] = 0.970861;
  phiscale[67][1] = 0.956949;
  phiscale[68][1] = 0.940985;
  phiscale[69][1] = 0.963053;
  phiscale[70][1] = 0.965262;
  phiscale[71][1] = 0.976937;
  phiscale[72][1] = 1.02766;
  phiscale[73][1] = 1.01104;
  phiscale[74][1] = 0.961209;
  phiscale[75][1] = 0.982048;
  phiscale[76][1] = 0.937176;
  phiscale[77][1] = 0.963477;
  phiscale[78][1] = 0.903683;
  phiscale[79][1] = 0.962511;
  phiscale[80][1] = 0.936425;
  phiscale[81][1] = 0.961669;
  phiscale[82][1] = 0.987199;
  phiscale[83][1] = 0.949005;
  phiscale[84][1] = 0.922117;
  phiscale[85][1] = 0.950696;
  phiscale[86][1] = 0.984253;
  phiscale[87][1] = 0.954781;
  phiscale[88][1] = 0.964978;
  phiscale[89][1] = 0.928218;
  phiscale[0][2] = 0.962744;
  phiscale[1][2] = 0.968236;
  phiscale[2][2] = 0.987835;
  phiscale[3][2] = 0.954738;
  phiscale[4][2] = 0.963859;
  phiscale[5][2] = 0.951937;
  phiscale[6][2] = 0.963344;
  phiscale[7][2] = 0.962387;
  phiscale[8][2] = 0.94448;
  phiscale[9][2] = 0.930095;
  phiscale[10][2] = 0.909768;
  phiscale[11][2] = 0.944394;
  phiscale[12][2] = 0.957912;
  phiscale[13][2] = 0.926358;
  phiscale[14][2] = 0.938448;
  phiscale[15][2] = 0.959608;
  phiscale[16][2] = 0.929638;
  phiscale[17][2] = 0.956948;
  phiscale[18][2] = 0.933224;
  phiscale[19][2] = 0.935572;
  phiscale[20][2] = 0.929563;
  phiscale[21][2] = 0.905672;
  phiscale[22][2] = 0.919862;
  phiscale[23][2] = 0.916125;
  phiscale[24][2] = 0.920462;
  phiscale[25][2] = 0.875855;
  phiscale[26][2] = 0.883922;
  phiscale[27][2] = 0.890918;
  phiscale[28][2] = 0.89703;
  phiscale[29][2] = 0.907915;
  phiscale[30][2] = 0.933024;
  phiscale[31][2] = 0.899936;
  phiscale[32][2] = 0.913392;
  phiscale[33][2] = 0.90232;
  phiscale[34][2] = 0.910455;
  phiscale[35][2] = 0.902662;
  phiscale[36][2] = 0.926537;
  phiscale[37][2] = 0.936633;
  phiscale[38][2] = 0.943717;
  phiscale[39][2] = 0.938017;
  phiscale[40][2] = 0.939619;
  phiscale[41][2] = 0.926241;
  phiscale[42][2] = 0.908418;
  phiscale[43][2] = 0.920894;
  phiscale[44][2] = 0.927385;
  phiscale[45][2] = 0.946894;
  phiscale[46][2] = 0.924502;
  phiscale[47][2] = 0.944671;
  phiscale[48][2] = 0.923533;
  phiscale[49][2] = 0.903496;
  phiscale[50][2] = 0.956848;
  phiscale[51][2] = 0.959836;
  phiscale[52][2] = 0.930753;
  phiscale[53][2] = 0.930292;
  phiscale[54][2] = 0.962006;
  phiscale[55][2] = 0.943203;
  phiscale[56][2] = 0.92875;
  phiscale[57][2] = 0.912189;
  phiscale[58][2] = 0.941786;
  phiscale[59][2] = 0.967458;
  phiscale[60][2] = 0.952121;
  phiscale[61][2] = 0.948033;
  phiscale[62][2] = 0.951697;
  phiscale[63][2] = 0.912986;
  phiscale[64][2] = 0.916005;
  phiscale[65][2] = 0.948351;
  phiscale[66][2] = 0.968604;
  phiscale[67][2] = 0.975172;
  phiscale[68][2] = 0.938047;
  phiscale[69][2] = 0.956017;
  phiscale[70][2] = 0.958354;
  phiscale[71][2] = 0.946423;
  phiscale[72][2] = 0.990133;
  phiscale[73][2] = 0.991984;
  phiscale[74][2] = 0.990906;
  phiscale[75][2] = 0.983481;
  phiscale[76][2] = 0.988379;
  phiscale[77][2] = 0.966493;
  phiscale[78][2] = 0.921903;
  phiscale[79][2] = 0.94546;
  phiscale[80][2] = 0.935957;
  phiscale[81][2] = 0.952414;
  phiscale[82][2] = 0.983684;
  phiscale[83][2] = 0.9571;
  phiscale[84][2] = 0.947858;
  phiscale[85][2] = 0.966653;
  phiscale[86][2] = 0.966577;
  phiscale[87][2] = 0.952644;
  phiscale[88][2] = 0.940709;
  phiscale[89][2] = 0.963862;
  phiscale[0][3] = 0.987794;
  phiscale[1][3] = 0.987147;
  phiscale[2][3] = 0.926399;
  phiscale[3][3] = 0.961158;
  phiscale[4][3] = 0.971941;
  phiscale[5][3] = 0.969951;
  phiscale[6][3] = 0.967806;
  phiscale[7][3] = 0.978516;
  phiscale[8][3] = 0.95474;
  phiscale[9][3] = 0.929099;
  phiscale[10][3] = 0.948314;
  phiscale[11][3] = 0.972973;
  phiscale[12][3] = 0.959863;
  phiscale[13][3] = 0.940294;
  phiscale[14][3] = 0.940651;
  phiscale[15][3] = 0.958237;
  phiscale[16][3] = 0.936998;
  phiscale[17][3] = 0.952006;
  phiscale[18][3] = 0.931194;
  phiscale[19][3] = 0.935452;
  phiscale[20][3] = 0.956579;
  phiscale[21][3] = 0.940007;
  phiscale[22][3] = 0.928022;
  phiscale[23][3] = 0.926449;
  phiscale[24][3] = 0.919303;
  phiscale[25][3] = 0.893119;
  phiscale[26][3] = 0.897961;
  phiscale[27][3] = 0.892992;
  phiscale[28][3] = 0.898598;
  phiscale[29][3] = 0.9214;
  phiscale[30][3] = 0.935744;
  phiscale[31][3] = 0.919061;
  phiscale[32][3] = 0.933211;
  phiscale[33][3] = 0.903891;
  phiscale[34][3] = 0.91554;
  phiscale[35][3] = 0.914643;
  phiscale[36][3] = 0.909571;
  phiscale[37][3] = 0.951164;
  phiscale[38][3] = 0.945412;
  phiscale[39][3] = 0.954608;
  phiscale[40][3] = 0.958747;
  phiscale[41][3] = 0.942217;
  phiscale[42][3] = 0.900071;
  phiscale[43][3] = 0.942896;
  phiscale[44][3] = 0.957917;
  phiscale[45][3] = 0.976359;
  phiscale[46][3] = 0.944075;
  phiscale[47][3] = 0.96731;
  phiscale[48][3] = 0.948096;
  phiscale[49][3] = 0.908108;
  phiscale[50][3] = 0.969574;
  phiscale[51][3] = 0.932715;
  phiscale[52][3] = 0.957889;
  phiscale[53][3] = 0.965732;
  phiscale[54][3] = 0.981088;
  phiscale[55][3] = 0.995867;
  phiscale[56][3] = 0.964963;
  phiscale[57][3] = 0.959919;
  phiscale[58][3] = 0.973134;
  phiscale[59][3] = 0.978695;
  phiscale[60][3] = 0.961887;
  phiscale[61][3] = 0.961397;
  phiscale[62][3] = 0.984122;
  phiscale[63][3] = 0.966031;
  phiscale[64][3] = 0.962322;
  phiscale[65][3] = 0.965317;
  phiscale[66][3] = 1.00935;
  phiscale[67][3] = 0.985189;
  phiscale[68][3] = 0.960067;
  phiscale[69][3] = 0.973611;
  phiscale[70][3] = 0.972342;
  phiscale[71][3] = 0.995384;
  phiscale[72][3] = 1.0381;
  phiscale[73][3] = 1.03394;
  phiscale[74][3] = 1.00726;
  phiscale[75][3] = 1.00078;
  phiscale[76][3] = 0.978509;
  phiscale[77][3] = 0.964068;
  phiscale[78][3] = 0.963376;
  phiscale[79][3] = 0.965202;
  phiscale[80][3] = 0.961356;
  phiscale[81][3] = 0.967696;
  phiscale[82][3] = 0.983539;
  phiscale[83][3] = 0.968968;
  phiscale[84][3] = 0.967558;
  phiscale[85][3] = 0.97386;
  phiscale[86][3] = 0.991842;
  phiscale[87][3] = 0.968669;
  phiscale[88][3] = 0.957027;
  phiscale[89][3] = 0.979929;
  phiscale[0][4] = 1.03479;
  phiscale[1][4] = 1.03132;
  phiscale[2][4] = 1.01542;
  phiscale[3][4] = 0.998314;
  phiscale[4][4] = 1.02657;
  phiscale[5][4] = 1.02967;
  phiscale[6][4] = 0.992283;
  phiscale[7][4] = 1.00198;
  phiscale[8][4] = 0.993084;
  phiscale[9][4] = 0.981729;
  phiscale[10][4] = 0.974627;
  phiscale[11][4] = 1.00631;
  phiscale[12][4] = 1.00988;
  phiscale[13][4] = 0.973337;
  phiscale[14][4] = 0.983055;
  phiscale[15][4] = 0.997373;
  phiscale[16][4] = 0.99266;
  phiscale[17][4] = 1.00054;
  phiscale[18][4] = 0.970996;
  phiscale[19][4] = 0.972716;
  phiscale[20][4] = 0.988278;
  phiscale[21][4] = 0.975554;
  phiscale[22][4] = 0.966111;
  phiscale[23][4] = 0.949356;
  phiscale[24][4] = 0.958193;
  phiscale[25][4] = 0.933206;
  phiscale[26][4] = 0.925199;
  phiscale[27][4] = 0.92394;
  phiscale[28][4] = 0.937992;
  phiscale[29][4] = 0.947246;
  phiscale[30][4] = 0.957649;
  phiscale[31][4] = 0.948422;
  phiscale[32][4] = 0.948078;
  phiscale[33][4] = 0.947662;
  phiscale[34][4] = 0.944964;
  phiscale[35][4] = 0.938989;
  phiscale[36][4] = 0.964;
  phiscale[37][4] = 0.970471;
  phiscale[38][4] = 0.978393;
  phiscale[39][4] = 0.978171;
  phiscale[40][4] = 0.993568;
  phiscale[41][4] = 0.967414;
  phiscale[42][4] = 0.968406;
  phiscale[43][4] = 0.964025;
  phiscale[44][4] = 0.975969;
  phiscale[45][4] = 1.00292;
  phiscale[46][4] = 0.975306;
  phiscale[47][4] = 1.00263;
  phiscale[48][4] = 0.988438;
  phiscale[49][4] = 0.990133;
  phiscale[50][4] = 1.00889;
  phiscale[51][4] = 0.990981;
  phiscale[52][4] = 0.976207;
  phiscale[53][4] = 0.965798;
  phiscale[54][4] = 1.01522;
  phiscale[55][4] = 1.0133;
  phiscale[56][4] = 1.02513;
  phiscale[57][4] = 1.02705;
  phiscale[58][4] = 1.00009;
  phiscale[59][4] = 1.02098;
  phiscale[60][4] = 1.0051;
  phiscale[61][4] = 0.987233;
  phiscale[62][4] = 1.0271;
  phiscale[63][4] = 0.997755;
  phiscale[64][4] = 0.999692;
  phiscale[65][4] = 1.00428;
  phiscale[66][4] = 1.03572;
  phiscale[67][4] = 1.0233;
  phiscale[68][4] = 1.01028;
  phiscale[69][4] = 1.00859;
  phiscale[70][4] = 1.01064;
  phiscale[71][4] = 1.04185;
  phiscale[72][4] = 1.07233;
  phiscale[73][4] = 1.07773;
  phiscale[74][4] = 1.05368;
  phiscale[75][4] = 1.03807;
  phiscale[76][4] = 1.01813;
  phiscale[77][4] = 1.00684;
  phiscale[78][4] = 0.980259;
  phiscale[79][4] = 1.00157;
  phiscale[80][4] = 0.998194;
  phiscale[81][4] = 0.999355;
  phiscale[82][4] = 1.0262;
  phiscale[83][4] = 0.994544;
  phiscale[84][4] = 0.995491;
  phiscale[85][4] = 1.0019;
  phiscale[86][4] = 1.03767;
  phiscale[87][4] = 1.01081;
  phiscale[88][4] = 0.986517;
  phiscale[89][4] = 1.01317;
  phiscale[0][5] = 1.00793;
  phiscale[1][5] = 1.0146;
  phiscale[2][5] = 0.983486;
  phiscale[3][5] = 0.988886;
  phiscale[4][5] = 1.01324;
  phiscale[5][5] = 1.00735;
  phiscale[6][5] = 0.991113;
  phiscale[7][5] = 1.00968;
  phiscale[8][5] = 0.990717;
  phiscale[9][5] = 0.976264;
  phiscale[10][5] = 0.975157;
  phiscale[11][5] = 0.996761;
  phiscale[12][5] = 1.01459;
  phiscale[13][5] = 0.983219;
  phiscale[14][5] = 0.994977;
  phiscale[15][5] = 0.991754;
  phiscale[16][5] = 0.978266;
  phiscale[17][5] = 0.997686;
  phiscale[18][5] = 0.96588;
  phiscale[19][5] = 0.976652;
  phiscale[20][5] = 0.988941;
  phiscale[21][5] = 0.960569;
  phiscale[22][5] = 0.960808;
  phiscale[23][5] = 0.938716;
  phiscale[24][5] = 0.940917;
  phiscale[25][5] = 0.916408;
  phiscale[26][5] = 0.911269;
  phiscale[27][5] = 0.907248;
  phiscale[28][5] = 0.920526;
  phiscale[29][5] = 0.931964;
  phiscale[30][5] = 0.953702;
  phiscale[31][5] = 0.928047;
  phiscale[32][5] = 0.929093;
  phiscale[33][5] = 0.933683;
  phiscale[34][5] = 0.932729;
  phiscale[35][5] = 0.931075;
  phiscale[36][5] = 0.943249;
  phiscale[37][5] = 0.962104;
  phiscale[38][5] = 0.957959;
  phiscale[39][5] = 0.966886;
  phiscale[40][5] = 0.969283;
  phiscale[41][5] = 0.956305;
  phiscale[42][5] = 0.942612;
  phiscale[43][5] = 0.955031;
  phiscale[44][5] = 0.965364;
  phiscale[45][5] = 0.986797;
  phiscale[46][5] = 0.967976;
  phiscale[47][5] = 0.978287;
  phiscale[48][5] = 0.96841;
  phiscale[49][5] = 0.976684;
  phiscale[50][5] = 1.00824;
  phiscale[51][5] = 0.989693;
  phiscale[52][5] = 0.994792;
  phiscale[53][5] = 0.999946;
  phiscale[54][5] = 1.00359;
  phiscale[55][5] = 1.01753;
  phiscale[56][5] = 1.02788;
  phiscale[57][5] = 1.00931;
  phiscale[58][5] = 1.03748;
  phiscale[59][5] = 1.02942;
  phiscale[60][5] = 0.996674;
  phiscale[61][5] = 1.01538;
  phiscale[62][5] = 1.02397;
  phiscale[63][5] = 1.01385;
  phiscale[64][5] = 1.024;
  phiscale[65][5] = 1.00126;
  phiscale[66][5] = 1.03374;
  phiscale[67][5] = 1.03196;
  phiscale[68][5] = 1.02336;
  phiscale[69][5] = 1.03348;
  phiscale[70][5] = 1.03293;
  phiscale[71][5] = 1.04117;
  phiscale[72][5] = 1.07982;
  phiscale[73][5] = 1.08375;
  phiscale[74][5] = 1.05868;
  phiscale[75][5] = 1.04345;
  phiscale[76][5] = 1.03456;
  phiscale[77][5] = 1.02341;
  phiscale[78][5] = 1.00026;
  phiscale[79][5] = 1.00516;
  phiscale[80][5] = 0.992589;
  phiscale[81][5] = 0.996006;
  phiscale[82][5] = 1.02293;
  phiscale[83][5] = 0.990703;
  phiscale[84][5] = 0.986053;
  phiscale[85][5] = 0.993369;
  phiscale[86][5] = 1.00889;
  phiscale[87][5] = 1.00669;
  phiscale[88][5] = 0.987497;
  phiscale[89][5] = 1.00253;
  phiscale[0][6] = 1.02731;
  phiscale[1][6] = 1.0406;
  phiscale[2][6] = 0.994241;
  phiscale[3][6] = 1.02427;
  phiscale[4][6] = 1.03308;
  phiscale[5][6] = 1.04346;
  phiscale[6][6] = 1.00812;
  phiscale[7][6] = 1.03361;
  phiscale[8][6] = 1.00317;
  phiscale[9][6] = 0.997491;
  phiscale[10][6] = 1.02213;
  phiscale[11][6] = 1.01924;
  phiscale[12][6] = 1.02554;
  phiscale[13][6] = 0.991703;
  phiscale[14][6] = 1.00515;
  phiscale[15][6] = 1.01845;
  phiscale[16][6] = 0.996485;
  phiscale[17][6] = 1.01683;
  phiscale[18][6] = 0.986422;
  phiscale[19][6] = 0.99701;
  phiscale[20][6] = 1.00656;
  phiscale[21][6] = 0.975031;
  phiscale[22][6] = 0.984748;
  phiscale[23][6] = 0.958899;
  phiscale[24][6] = 0.95318;
  phiscale[25][6] = 0.945882;
  phiscale[26][6] = 0.937938;
  phiscale[27][6] = 0.92109;
  phiscale[28][6] = 0.929799;
  phiscale[29][6] = 0.944901;
  phiscale[30][6] = 0.973326;
  phiscale[31][6] = 0.947088;
  phiscale[32][6] = 0.950689;
  phiscale[33][6] = 0.949488;
  phiscale[34][6] = 0.945845;
  phiscale[35][6] = 0.938223;
  phiscale[36][6] = 0.956736;
  phiscale[37][6] = 0.967982;
  phiscale[38][6] = 0.97853;
  phiscale[39][6] = 0.976224;
  phiscale[40][6] = 0.986344;
  phiscale[41][6] = 0.968759;
  phiscale[42][6] = 0.960719;
  phiscale[43][6] = 0.979602;
  phiscale[44][6] = 0.973673;
  phiscale[45][6] = 1.0053;
  phiscale[46][6] = 0.967541;
  phiscale[47][6] = 1.00028;
  phiscale[48][6] = 0.992634;
  phiscale[49][6] = 0.965707;
  phiscale[50][6] = 1.01571;
  phiscale[51][6] = 1.00435;
  phiscale[52][6] = 1.02529;
  phiscale[53][6] = 1.00509;
  phiscale[54][6] = 1.04223;
  phiscale[55][6] = 1.05097;
  phiscale[56][6] = 1.039;
  phiscale[57][6] = 1.02508;
  phiscale[58][6] = 1.06977;
  phiscale[59][6] = 1.05493;
  phiscale[60][6] = 1.03653;
  phiscale[61][6] = 1.04124;
  phiscale[62][6] = 1.05847;
  phiscale[63][6] = 1.05553;
  phiscale[64][6] = 1.04456;
  phiscale[65][6] = 1.06792;
  phiscale[66][6] = 1.09019;
  phiscale[67][6] = 1.07576;
  phiscale[68][6] = 1.06413;
  phiscale[69][6] = 1.06683;
  phiscale[70][6] = 1.07147;
  phiscale[71][6] = 1.08243;
  phiscale[72][6] = 1.12099;
  phiscale[73][6] = 1.11348;
  phiscale[74][6] = 1.09752;
  phiscale[75][6] = 1.07734;
  phiscale[76][6] = 1.06748;
  phiscale[77][6] = 1.05381;
  phiscale[78][6] = 1.0353;
  phiscale[79][6] = 1.02817;
  phiscale[80][6] = 1.00546;
  phiscale[81][6] = 1.00717;
  phiscale[82][6] = 1.02958;
  phiscale[83][6] = 1.01589;
  phiscale[84][6] = 1.01348;
  phiscale[85][6] = 1.03639;
  phiscale[86][6] = 1.02613;
  phiscale[87][6] = 1.0099;
  phiscale[88][6] = 0.997602;
  phiscale[89][6] = 1.0081;
  phiscale[0][7] = 1.02825;
  phiscale[1][7] = 1.03708;
  phiscale[2][7] = 1.02834;
  phiscale[3][7] = 1.01697;
  phiscale[4][7] = 1.04575;
  phiscale[5][7] = 1.04302;
  phiscale[6][7] = 1.02678;
  phiscale[7][7] = 1.03781;
  phiscale[8][7] = 1.01394;
  phiscale[9][7] = 0.997676;
  phiscale[10][7] = 1.01142;
  phiscale[11][7] = 1.03645;
  phiscale[12][7] = 1.03594;
  phiscale[13][7] = 1.01234;
  phiscale[14][7] = 1.02392;
  phiscale[15][7] = 1.02677;
  phiscale[16][7] = 1.01465;
  phiscale[17][7] = 1.03232;
  phiscale[18][7] = 1.00286;
  phiscale[19][7] = 0.997639;
  phiscale[20][7] = 1.01477;
  phiscale[21][7] = 0.99223;
  phiscale[22][7] = 0.985837;
  phiscale[23][7] = 0.964253;
  phiscale[24][7] = 0.95424;
  phiscale[25][7] = 0.93627;
  phiscale[26][7] = 0.932802;
  phiscale[27][7] = 0.924114;
  phiscale[28][7] = 0.925617;
  phiscale[29][7] = 0.941879;
  phiscale[30][7] = 0.967933;
  phiscale[31][7] = 0.938107;
  phiscale[32][7] = 0.953124;
  phiscale[33][7] = 0.941306;
  phiscale[34][7] = 0.947676;
  phiscale[35][7] = 0.933013;
  phiscale[36][7] = 0.946559;
  phiscale[37][7] = 0.968474;
  phiscale[38][7] = 0.972704;
  phiscale[39][7] = 0.976576;
  phiscale[40][7] = 0.986255;
  phiscale[41][7] = 0.974583;
  phiscale[42][7] = 0.955271;
  phiscale[43][7] = 0.977918;
  phiscale[44][7] = 0.990417;
  phiscale[45][7] = 0.994721;
  phiscale[46][7] = 0.983102;
  phiscale[47][7] = 0.996746;
  phiscale[48][7] = 0.971459;
  phiscale[49][7] = 0.935865;
  phiscale[50][7] = 1.03571;
  phiscale[51][7] = 1.03211;
  phiscale[52][7] = 1.01899;
  phiscale[53][7] = 1.01937;
  phiscale[54][7] = 1.03731;
  phiscale[55][7] = 1.06579;
  phiscale[56][7] = 1.05733;
  phiscale[57][7] = 1.03967;
  phiscale[58][7] = 1.06825;
  phiscale[59][7] = 1.08051;
  phiscale[60][7] = 1.05746;
  phiscale[61][7] = 1.05815;
  phiscale[62][7] = 1.08345;
  phiscale[63][7] = 1.07144;
  phiscale[64][7] = 1.05635;
  phiscale[65][7] = 1.06699;
  phiscale[66][7] = 1.09864;
  phiscale[67][7] = 1.08891;
  phiscale[68][7] = 1.08466;
  phiscale[69][7] = 1.08825;
  phiscale[70][7] = 1.08409;
  phiscale[71][7] = 1.1061;
  phiscale[72][7] = 1.13566;
  phiscale[73][7] = 1.10155;
  phiscale[74][7] = 1.10043;
  phiscale[75][7] = 1.08027;
  phiscale[76][7] = 1.06131;
  phiscale[77][7] = 1.04556;
  phiscale[78][7] = 1.02172;
  phiscale[79][7] = 1.04103;
  phiscale[80][7] = 1.0161;
  phiscale[81][7] = 1.01957;
  phiscale[82][7] = 1.03213;
  phiscale[83][7] = 1.01323;
  phiscale[84][7] = 1.00874;
  phiscale[85][7] = 1.01865;
  phiscale[86][7] = 1.03355;
  phiscale[87][7] = 1.00922;
  phiscale[88][7] = 0.997355;
  phiscale[89][7] = 1.01721;
  phiscale[0][8] = 1.0152;
  phiscale[1][8] = 1.02942;
  phiscale[2][8] = 0.985572;
  phiscale[3][8] = 1.01253;
  phiscale[4][8] = 1.02215;
  phiscale[5][8] = 1.04269;
  phiscale[6][8] = 1.01508;
  phiscale[7][8] = 1.03454;
  phiscale[8][8] = 1.00203;
  phiscale[9][8] = 0.995241;
  phiscale[10][8] = 1.00939;
  phiscale[11][8] = 1.01173;
  phiscale[12][8] = 1.00769;
  phiscale[13][8] = 0.999926;
  phiscale[14][8] = 1.01162;
  phiscale[15][8] = 1.01751;
  phiscale[16][8] = 0.999405;
  phiscale[17][8] = 1.01172;
  phiscale[18][8] = 0.992602;
  phiscale[19][8] = 0.995862;
  phiscale[20][8] = 0.999937;
  phiscale[21][8] = 0.982535;
  phiscale[22][8] = 0.968637;
  phiscale[23][8] = 0.94445;
  phiscale[24][8] = 0.940227;
  phiscale[25][8] = 0.913647;
  phiscale[26][8] = 0.915205;
  phiscale[27][8] = 0.89235;
  phiscale[28][8] = 0.91458;
  phiscale[29][8] = 0.925375;
  phiscale[30][8] = 0.945705;
  phiscale[31][8] = 0.934649;
  phiscale[32][8] = 0.937861;
  phiscale[33][8] = 0.933295;
  phiscale[34][8] = 0.917147;
  phiscale[35][8] = 0.922435;
  phiscale[36][8] = 0.926805;
  phiscale[37][8] = 0.950925;
  phiscale[38][8] = 0.962165;
  phiscale[39][8] = 0.963178;
  phiscale[40][8] = 0.981573;
  phiscale[41][8] = 0.958654;
  phiscale[42][8] = 0.97082;
  phiscale[43][8] = 0.959631;
  phiscale[44][8] = 0.963387;
  phiscale[45][8] = 0.988221;
  phiscale[46][8] = 0.97261;
  phiscale[47][8] = 0.980705;
  phiscale[48][8] = 0.977725;
  phiscale[49][8] = 0.981512;
  phiscale[50][8] = 1.01109;
  phiscale[51][8] = 1.00432;
  phiscale[52][8] = 1.01343;
  phiscale[53][8] = 1.01408;
  phiscale[54][8] = 1.03918;
  phiscale[55][8] = 1.05154;
  phiscale[56][8] = 1.04364;
  phiscale[57][8] = 1.06335;
  phiscale[58][8] = 1.06097;
  phiscale[59][8] = 1.07421;
  phiscale[60][8] = 1.06173;
  phiscale[61][8] = 1.05687;
  phiscale[62][8] = 1.07908;
  phiscale[63][8] = 1.0642;
  phiscale[64][8] = 1.05737;
  phiscale[65][8] = 1.0819;
  phiscale[66][8] = 1.11396;
  phiscale[67][8] = 1.09762;
  phiscale[68][8] = 1.07319;
  phiscale[69][8] = 1.10362;
  phiscale[70][8] = 1.08647;
  phiscale[71][8] = 1.09291;
  phiscale[72][8] = 1.12696;
  phiscale[73][8] = 1.10733;
  phiscale[74][8] = 1.09018;
  phiscale[75][8] = 1.0599;
  phiscale[76][8] = 1.05004;
  phiscale[77][8] = 1.03903;
  phiscale[78][8] = 1.02628;
  phiscale[79][8] = 1.03628;
  phiscale[80][8] = 1.01461;
  phiscale[81][8] = 1.00736;
  phiscale[82][8] = 1.02583;
  phiscale[83][8] = 0.989613;
  phiscale[84][8] = 0.993896;
  phiscale[85][8] = 1.00182;
  phiscale[86][8] = 1.01037;
  phiscale[87][8] = 1.00198;
  phiscale[88][8] = 0.990354;
  phiscale[89][8] = 1.00193;
  phiscale[0][9] = 1.06127;
  phiscale[1][9] = 1.09003;
  phiscale[2][9] = 1.07286;
  phiscale[3][9] = 1.07831;
  phiscale[4][9] = 1.09432;
  phiscale[5][9] = 1.09833;
  phiscale[6][9] = 1.08445;
  phiscale[7][9] = 1.09625;
  phiscale[8][9] = 1.07188;
  phiscale[9][9] = 1.0631;
  phiscale[10][9] = 1.05567;
  phiscale[11][9] = 1.08234;
  phiscale[12][9] = 1.0851;
  phiscale[13][9] = 1.05887;
  phiscale[14][9] = 1.07198;
  phiscale[15][9] = 1.07402;
  phiscale[16][9] = 1.06566;
  phiscale[17][9] = 1.07314;
  phiscale[18][9] = 1.04131;
  phiscale[19][9] = 1.05322;
  phiscale[20][9] = 1.06025;
  phiscale[21][9] = 1.0511;
  phiscale[22][9] = 1.02637;
  phiscale[23][9] = 1.00702;
  phiscale[24][9] = 0.995874;
  phiscale[25][9] = 0.961624;
  phiscale[26][9] = 0.960818;
  phiscale[27][9] = 0.94267;
  phiscale[28][9] = 0.94804;
  phiscale[29][9] = 0.965843;
  phiscale[30][9] = 0.998068;
  phiscale[31][9] = 0.96609;
  phiscale[32][9] = 0.980801;
  phiscale[33][9] = 0.962131;
  phiscale[34][9] = 0.968221;
  phiscale[35][9] = 0.96725;
  phiscale[36][9] = 0.974776;
  phiscale[37][9] = 1.00226;
  phiscale[38][9] = 0.998606;
  phiscale[39][9] = 1.00801;
  phiscale[40][9] = 1.01182;
  phiscale[41][9] = 0.993283;
  phiscale[42][9] = 0.98955;
  phiscale[43][9] = 1.00372;
  phiscale[44][9] = 1.01267;
  phiscale[45][9] = 1.02436;
  phiscale[46][9] = 1.01523;
  phiscale[47][9] = 1.03405;
  phiscale[48][9] = 1.0343;
  phiscale[49][9] = 0.993256;
  phiscale[50][9] = 1.0857;
  phiscale[51][9] = 1.06796;
  phiscale[52][9] = 1.07578;
  phiscale[53][9] = 1.06764;
  phiscale[54][9] = 1.1043;
  phiscale[55][9] = 1.12782;
  phiscale[56][9] = 1.11739;
  phiscale[57][9] = 1.10687;
  phiscale[58][9] = 1.13576;
  phiscale[59][9] = 1.13361;
  phiscale[60][9] = 1.11387;
  phiscale[61][9] = 1.12756;
  phiscale[62][9] = 1.14286;
  phiscale[63][9] = 1.13035;
  phiscale[64][9] = 1.13981;
  phiscale[65][9] = 1.14405;
  phiscale[66][9] = 1.17879;
  phiscale[67][9] = 1.13626;
  phiscale[68][9] = 1.12439;
  phiscale[69][9] = 1.13084;
  phiscale[70][9] = 1.1164;
  phiscale[71][9] = 1.14558;
  phiscale[72][9] = 1.15872;
  phiscale[73][9] = 1.12942;
  phiscale[74][9] = 1.1282;
  phiscale[75][9] = 1.09076;
  phiscale[76][9] = 1.10039;
  phiscale[77][9] = 1.06331;
  phiscale[78][9] = 1.0409;
  phiscale[79][9] = 1.04297;
  phiscale[80][9] = 1.03209;
  phiscale[81][9] = 1.02696;
  phiscale[82][9] = 1.03767;
  phiscale[83][9] = 1.01207;
  phiscale[84][9] = 1.00439;
  phiscale[85][9] = 1.01389;
  phiscale[86][9] = 1.03773;
  phiscale[87][9] = 1.01711;
  phiscale[88][9] = 1.01801;
  phiscale[89][9] = 1.04554;
  phiscale[0][10] = 1.07205;
  phiscale[1][10] = 1.11765;
  phiscale[2][10] = 1.0686;
  phiscale[3][10] = 1.10075;
  phiscale[4][10] = 1.1171;
  phiscale[5][10] = 1.10962;
  phiscale[6][10] = 1.09173;
  phiscale[7][10] = 1.11042;
  phiscale[8][10] = 1.12456;
  phiscale[9][10] = 1.06544;
  phiscale[10][10] = 1.07606;
  phiscale[11][10] = 1.10121;
  phiscale[12][10] = 1.10567;
  phiscale[13][10] = 1.07496;
  phiscale[14][10] = 1.08953;
  phiscale[15][10] = 1.08052;
  phiscale[16][10] = 1.08789;
  phiscale[17][10] = 1.08583;
  phiscale[18][10] = 1.07895;
  phiscale[19][10] = 1.08765;
  phiscale[20][10] = 1.08382;
  phiscale[21][10] = 1.05601;
  phiscale[22][10] = 1.03913;
  phiscale[23][10] = 1.01212;
  phiscale[24][10] = 1.00518;
  phiscale[25][10] = 0.965911;
  phiscale[26][10] = 0.95409;
  phiscale[27][10] = 0.933396;
  phiscale[28][10] = 0.957111;
  phiscale[29][10] = 0.961018;
  phiscale[30][10] = 0.998547;
  phiscale[31][10] = 0.964817;
  phiscale[32][10] = 0.97026;
  phiscale[33][10] = 0.954038;
  phiscale[34][10] = 0.974994;
  phiscale[35][10] = 0.965788;
  phiscale[36][10] = 0.985841;
  phiscale[37][10] = 0.989116;
  phiscale[38][10] = 0.994092;
  phiscale[39][10] = 0.992494;
  phiscale[40][10] = 1.03559;
  phiscale[41][10] = 1.01504;
  phiscale[42][10] = 0.959925;
  phiscale[43][10] = 0.985876;
  phiscale[44][10] = 1.02605;
  phiscale[45][10] = 1.0231;
  phiscale[46][10] = 1.04125;
  phiscale[47][10] = 1.05077;
  phiscale[48][10] = 1.01309;
  phiscale[49][10] = 1.00869;
  phiscale[50][10] = 1.08454;
  phiscale[51][10] = 1.0744;
  phiscale[52][10] = 1.08205;
  phiscale[53][10] = 1.08699;
  phiscale[54][10] = 1.09914;
  phiscale[55][10] = 1.12222;
  phiscale[56][10] = 1.1482;
  phiscale[57][10] = 1.17036;
  phiscale[58][10] = 1.16502;
  phiscale[59][10] = 1.16622;
  phiscale[60][10] = 1.1455;
  phiscale[61][10] = 1.1545;
  phiscale[62][10] = 1.16757;
  phiscale[63][10] = 1.15305;
  phiscale[64][10] = 1.13882;
  phiscale[65][10] = 1.17146;
  phiscale[66][10] = 1.19605;
  phiscale[67][10] = 1.13314;
  phiscale[68][10] = 1.1606;
  phiscale[69][10] = 1.14638;
  phiscale[70][10] = 1.13672;
  phiscale[71][10] = 1.15801;
  phiscale[72][10] = 1.16374;
  phiscale[73][10] = 1.14097;
  phiscale[74][10] = 1.14194;
  phiscale[75][10] = 1.12196;
  phiscale[76][10] = 1.10052;
  phiscale[77][10] = 1.09457;
  phiscale[78][10] = 1.05275;
  phiscale[79][10] = 1.0404;
  phiscale[80][10] = 1.04067;
  phiscale[81][10] = 1.0279;
  phiscale[82][10] = 1.05961;
  phiscale[83][10] = 1.02525;
  phiscale[84][10] = 1.00882;
  phiscale[85][10] = 1.02646;
  phiscale[86][10] = 1.048;
  phiscale[87][10] = 1.02508;
  phiscale[88][10] = 1.01068;
  phiscale[89][10] = 1.06141;

  Int_t runIndex = RunIndex(run);
  Int_t phiIndex = PhiIndex(phi);

  return phiscale[phiIndex][runIndex];
}

//Phi Correction (negative charge)
//V3.0 02/24/2004 SG 
Double_t PhiCorrNeg(Int_t run, Double_t phi) {

  // corr = phiscale(index phi, index run)

  Double_t phiscale[90][11]; //90 bins in phi * 11 run ranges

  phiscale[0][0] = 0.998376;
  phiscale[1][0] = 0.981414;
  phiscale[2][0] = 0.971507;
  phiscale[3][0] = 0.921157;
  phiscale[4][0] = 0.985675;
  phiscale[5][0] = 0.935106;
  phiscale[6][0] = 1.00578;
  phiscale[7][0] = 0.939247;
  phiscale[8][0] = 0.964935;
  phiscale[9][0] = 0.992724;
  phiscale[10][0] = 0.973312;
  phiscale[11][0] = 0.949487;
  phiscale[12][0] = 0.947539;
  phiscale[13][0] = 0.991065;
  phiscale[14][0] = 0.956782;
  phiscale[15][0] = 0.957009;
  phiscale[16][0] = 0.972543;
  phiscale[17][0] = 0.938433;
  phiscale[18][0] = 0.971758;
  phiscale[19][0] = 0.973065;
  phiscale[20][0] = 0.94632;
  phiscale[21][0] = 0.961519;
  phiscale[22][0] = 0.959597;
  phiscale[23][0] = 0.972474;
  phiscale[24][0] = 0.923774;
  phiscale[25][0] = 0.921746;
  phiscale[26][0] = 0.935283;
  phiscale[27][0] = 0.928241;
  phiscale[28][0] = 0.915408;
  phiscale[29][0] = 0.933155;
  phiscale[30][0] = 0.935102;
  phiscale[31][0] = 0.955183;
  phiscale[32][0] = 0.956867;
  phiscale[33][0] = 0.933709;
  phiscale[34][0] = 0.938234;
  phiscale[35][0] = 0.953714;
  phiscale[36][0] = 0.937717;
  phiscale[37][0] = 0.969798;
  phiscale[38][0] = 0.959352;
  phiscale[39][0] = 0.970133;
  phiscale[40][0] = 0.991218;
  phiscale[41][0] = 0.979924;
  phiscale[42][0] = 0.960888;
  phiscale[43][0] = 0.945052;
  phiscale[44][0] = 0.938818;
  phiscale[45][0] = 0.965865;
  phiscale[46][0] = 0.979375;
  phiscale[47][0] = 0.976497;
  phiscale[48][0] = 0.9865;
  phiscale[49][0] = 0.968807;
  phiscale[50][0] = 0.975379;
  phiscale[51][0] = 0.984383;
  phiscale[52][0] = 0.961926;
  phiscale[53][0] = 0.991054;
  phiscale[54][0] = 0.949938;
  phiscale[55][0] = 0.954334;
  phiscale[56][0] = 0.955987;
  phiscale[57][0] = 0.985484;
  phiscale[58][0] = 1.00671;
  phiscale[59][0] = 0.98074;
  phiscale[60][0] = 0.958437;
  phiscale[61][0] = 0.978307;
  phiscale[62][0] = 0.961412;
  phiscale[63][0] = 0.968028;
  phiscale[64][0] = 0.970833;
  phiscale[65][0] = 0.950994;
  phiscale[66][0] = 0.94324;
  phiscale[67][0] = 0.991175;
  phiscale[68][0] = 1.00712;
  phiscale[69][0] = 0.95685;
  phiscale[70][0] = 0.9822;
  phiscale[71][0] = 1.00186;
  phiscale[72][0] = 0.979935;
  phiscale[73][0] = 1.02498;
  phiscale[74][0] = 1.04945;
  phiscale[75][0] = 1.02076;
  phiscale[76][0] = 1.01467;
  phiscale[77][0] = 0.985253;
  phiscale[78][0] = 0.972927;
  phiscale[79][0] = 0.949129;
  phiscale[80][0] = 0.956824;
  phiscale[81][0] = 0.968626;
  phiscale[82][0] = 0.992901;
  phiscale[83][0] = 0.962966;
  phiscale[84][0] = 0.997452;
  phiscale[85][0] = 0.976229;
  phiscale[86][0] = 1.00792;
  phiscale[87][0] = 0.982813;
  phiscale[88][0] = 1.00191;
  phiscale[89][0] = 0.962759;
  phiscale[0][1] = 0.9442;
  phiscale[1][1] = 0.944782;
  phiscale[2][1] = 1.05388;
  phiscale[3][1] = 0.919985;
  phiscale[4][1] = 0.949282;
  phiscale[5][1] = 0.958502;
  phiscale[6][1] = 0.982345;
  phiscale[7][1] = 0.930565;
  phiscale[8][1] = 0.941385;
  phiscale[9][1] = 0.962732;
  phiscale[10][1] = 0.951816;
  phiscale[11][1] = 0.888298;
  phiscale[12][1] = 0.931806;
  phiscale[13][1] = 0.979049;
  phiscale[14][1] = 0.924174;
  phiscale[15][1] = 0.943276;
  phiscale[16][1] = 0.958807;
  phiscale[17][1] = 0.915802;
  phiscale[18][1] = 0.966296;
  phiscale[19][1] = 0.988799;
  phiscale[20][1] = 0.947501;
  phiscale[21][1] = 0.983206;
  phiscale[22][1] = 0.956257;
  phiscale[23][1] = 0.948475;
  phiscale[24][1] = 0.930598;
  phiscale[25][1] = 0.910531;
  phiscale[26][1] = 0.898929;
  phiscale[27][1] = 0.898958;
  phiscale[28][1] = 0.897777;
  phiscale[29][1] = 0.921664;
  phiscale[30][1] = 0.916759;
  phiscale[31][1] = 0.931823;
  phiscale[32][1] = 0.951151;
  phiscale[33][1] = 0.914622;
  phiscale[34][1] = 0.899382;
  phiscale[35][1] = 0.9134;
  phiscale[36][1] = 0.906709;
  phiscale[37][1] = 0.914244;
  phiscale[38][1] = 0.943761;
  phiscale[39][1] = 0.924433;
  phiscale[40][1] = 0.942176;
  phiscale[41][1] = 0.968607;
  phiscale[42][1] = 0.930591;
  phiscale[43][1] = 0.894358;
  phiscale[44][1] = 0.928436;
  phiscale[45][1] = 0.886051;
  phiscale[46][1] = 0.96685;
  phiscale[47][1] = 0.967146;
  phiscale[48][1] = 0.94611;
  phiscale[49][1] = 0.919834;
  phiscale[50][1] = 0.987088;
  phiscale[51][1] = 0.955131;
  phiscale[52][1] = 0.918156;
  phiscale[53][1] = 0.959506;
  phiscale[54][1] = 0.95657;
  phiscale[55][1] = 0.929944;
  phiscale[56][1] = 0.96211;
  phiscale[57][1] = 0.964038;
  phiscale[58][1] = 1.00675;
  phiscale[59][1] = 0.951446;
  phiscale[60][1] = 0.957159;
  phiscale[61][1] = 0.982671;
  phiscale[62][1] = 0.946238;
  phiscale[63][1] = 1.02279;
  phiscale[64][1] = 0.955571;
  phiscale[65][1] = 0.941376;
  phiscale[66][1] = 0.944266;
  phiscale[67][1] = 0.941517;
  phiscale[68][1] = 0.987761;
  phiscale[69][1] = 0.940974;
  phiscale[70][1] = 0.954025;
  phiscale[71][1] = 0.94775;
  phiscale[72][1] = 0.959033;
  phiscale[73][1] = 1.00857;
  phiscale[74][1] = 1.00098;
  phiscale[75][1] = 0.990933;
  phiscale[76][1] = 0.964247;
  phiscale[77][1] = 0.960272;
  phiscale[78][1] = 0.936804;
  phiscale[79][1] = 0.911931;
  phiscale[80][1] = 0.946468;
  phiscale[81][1] = 0.941211;
  phiscale[82][1] = 0.933433;
  phiscale[83][1] = 0.984777;
  phiscale[84][1] = 0.929661;
  phiscale[85][1] = 0.953882;
  phiscale[86][1] = 0.940275;
  phiscale[87][1] = 0.955159;
  phiscale[88][1] = 0.98892;
  phiscale[89][1] = 0.97434;
  phiscale[0][2] = 0.966418;
  phiscale[1][2] = 0.964245;
  phiscale[2][2] = 0.989605;
  phiscale[3][2] = 0.946878;
  phiscale[4][2] = 0.962564;
  phiscale[5][2] = 0.961493;
  phiscale[6][2] = 0.997318;
  phiscale[7][2] = 0.955573;
  phiscale[8][2] = 0.956623;
  phiscale[9][2] = 0.968692;
  phiscale[10][2] = 0.934843;
  phiscale[11][2] = 0.925676;
  phiscale[12][2] = 0.925644;
  phiscale[13][2] = 0.970753;
  phiscale[14][2] = 0.962225;
  phiscale[15][2] = 0.947116;
  phiscale[16][2] = 0.948138;
  phiscale[17][2] = 0.934178;
  phiscale[18][2] = 0.966974;
  phiscale[19][2] = 0.933544;
  phiscale[20][2] = 0.924313;
  phiscale[21][2] = 0.94167;
  phiscale[22][2] = 0.936146;
  phiscale[23][2] = 0.924686;
  phiscale[24][2] = 0.927072;
  phiscale[25][2] = 0.914487;
  phiscale[26][2] = 0.906063;
  phiscale[27][2] = 0.891303;
  phiscale[28][2] = 0.887186;
  phiscale[29][2] = 0.907739;
  phiscale[30][2] = 0.901968;
  phiscale[31][2] = 0.924889;
  phiscale[32][2] = 0.929866;
  phiscale[33][2] = 0.919656;
  phiscale[34][2] = 0.921306;
  phiscale[35][2] = 0.879987;
  phiscale[36][2] = 0.908257;
  phiscale[37][2] = 0.934795;
  phiscale[38][2] = 0.923844;
  phiscale[39][2] = 0.929282;
  phiscale[40][2] = 0.948216;
  phiscale[41][2] = 0.946911;
  phiscale[42][2] = 0.968098;
  phiscale[43][2] = 0.917672;
  phiscale[44][2] = 0.916995;
  phiscale[45][2] = 0.934753;
  phiscale[46][2] = 0.954667;
  phiscale[47][2] = 0.942242;
  phiscale[48][2] = 0.925083;
  phiscale[49][2] = 0.955953;
  phiscale[50][2] = 0.941901;
  phiscale[51][2] = 0.952057;
  phiscale[52][2] = 0.972069;
  phiscale[53][2] = 0.963396;
  phiscale[54][2] = 0.92169;
  phiscale[55][2] = 0.950559;
  phiscale[56][2] = 0.970784;
  phiscale[57][2] = 0.937461;
  phiscale[58][2] = 0.979782;
  phiscale[59][2] = 0.936586;
  phiscale[60][2] = 0.950578;
  phiscale[61][2] = 0.945083;
  phiscale[62][2] = 0.950885;
  phiscale[63][2] = 0.964905;
  phiscale[64][2] = 0.961398;
  phiscale[65][2] = 0.939885;
  phiscale[66][2] = 0.946904;
  phiscale[67][2] = 0.960158;
  phiscale[68][2] = 0.973714;
  phiscale[69][2] = 0.955901;
  phiscale[70][2] = 0.969602;
  phiscale[71][2] = 0.971469;
  phiscale[72][2] = 0.954681;
  phiscale[73][2] = 1.00066;
  phiscale[74][2] = 1.04256;
  phiscale[75][2] = 1.00084;
  phiscale[76][2] = 0.967803;
  phiscale[77][2] = 0.979525;
  phiscale[78][2] = 0.955522;
  phiscale[79][2] = 0.947087;
  phiscale[80][2] = 0.934297;
  phiscale[81][2] = 0.933536;
  phiscale[82][2] = 0.95837;
  phiscale[83][2] = 0.978413;
  phiscale[84][2] = 0.955851;
  phiscale[85][2] = 0.93966;
  phiscale[86][2] = 0.947012;
  phiscale[87][2] = 0.967499;
  phiscale[88][2] = 0.983927;
  phiscale[89][2] = 0.983277;
  phiscale[0][3] = 0.965005;
  phiscale[1][3] = 0.985449;
  phiscale[2][3] = 0.991375;
  phiscale[3][3] = 0.960516;
  phiscale[4][3] = 0.974327;
  phiscale[5][3] = 0.976429;
  phiscale[6][3] = 1.00244;
  phiscale[7][3] = 0.973423;
  phiscale[8][3] = 0.968253;
  phiscale[9][3] = 0.96653;
  phiscale[10][3] = 0.937257;
  phiscale[11][3] = 0.940596;
  phiscale[12][3] = 0.949652;
  phiscale[13][3] = 0.980453;
  phiscale[14][3] = 0.965471;
  phiscale[15][3] = 0.944693;
  phiscale[16][3] = 0.967965;
  phiscale[17][3] = 0.940574;
  phiscale[18][3] = 0.976954;
  phiscale[19][3] = 0.947866;
  phiscale[20][3] = 0.932218;
  phiscale[21][3] = 0.94259;
  phiscale[22][3] = 0.965046;
  phiscale[23][3] = 0.939968;
  phiscale[24][3] = 0.916412;
  phiscale[25][3] = 0.93332;
  phiscale[26][3] = 0.919988;
  phiscale[27][3] = 0.905102;
  phiscale[28][3] = 0.908464;
  phiscale[29][3] = 0.914644;
  phiscale[30][3] = 0.907568;
  phiscale[31][3] = 0.934449;
  phiscale[32][3] = 0.955457;
  phiscale[33][3] = 0.934077;
  phiscale[34][3] = 0.931215;
  phiscale[35][3] = 0.912752;
  phiscale[36][3] = 0.917747;
  phiscale[37][3] = 0.932262;
  phiscale[38][3] = 0.92949;
  phiscale[39][3] = 0.946847;
  phiscale[40][3] = 0.950908;
  phiscale[41][3] = 0.964686;
  phiscale[42][3] = 0.935517;
  phiscale[43][3] = 0.937816;
  phiscale[44][3] = 0.928775;
  phiscale[45][3] = 0.959223;
  phiscale[46][3] = 0.967295;
  phiscale[47][3] = 0.942993;
  phiscale[48][3] = 0.965636;
  phiscale[49][3] = 0.949432;
  phiscale[50][3] = 0.960852;
  phiscale[51][3] = 0.967078;
  phiscale[52][3] = 0.97588;
  phiscale[53][3] = 0.953466;
  phiscale[54][3] = 0.967758;
  phiscale[55][3] = 0.967616;
  phiscale[56][3] = 0.985651;
  phiscale[57][3] = 0.966892;
  phiscale[58][3] = 1.00292;
  phiscale[59][3] = 0.973899;
  phiscale[60][3] = 0.986005;
  phiscale[61][3] = 0.977844;
  phiscale[62][3] = 0.974313;
  phiscale[63][3] = 0.979924;
  phiscale[64][3] = 0.97139;
  phiscale[65][3] = 0.949296;
  phiscale[66][3] = 0.959119;
  phiscale[67][3] = 0.980217;
  phiscale[68][3] = 0.979498;
  phiscale[69][3] = 0.97315;
  phiscale[70][3] = 0.978272;
  phiscale[71][3] = 1.0067;
  phiscale[72][3] = 0.989686;
  phiscale[73][3] = 1.0331;
  phiscale[74][3] = 1.04254;
  phiscale[75][3] = 1.01531;
  phiscale[76][3] = 1.00132;
  phiscale[77][3] = 0.994371;
  phiscale[78][3] = 0.97643;
  phiscale[79][3] = 0.955948;
  phiscale[80][3] = 0.960772;
  phiscale[81][3] = 0.971661;
  phiscale[82][3] = 0.988556;
  phiscale[83][3] = 0.978384;
  phiscale[84][3] = 0.975106;
  phiscale[85][3] = 0.94791;
  phiscale[86][3] = 0.981169;
  phiscale[87][3] = 0.984478;
  phiscale[88][3] = 1.00783;
  phiscale[89][3] = 0.972637;
  phiscale[0][4] = 0.999869;
  phiscale[1][4] = 1.01584;
  phiscale[2][4] = 1.0253;
  phiscale[3][4] = 0.9926;
  phiscale[4][4] = 1.00428;
  phiscale[5][4] = 1.00862;
  phiscale[6][4] = 1.03247;
  phiscale[7][4] = 1.01803;
  phiscale[8][4] = 1.00666;
  phiscale[9][4] = 0.999054;
  phiscale[10][4] = 0.972024;
  phiscale[11][4] = 0.985414;
  phiscale[12][4] = 0.992628;
  phiscale[13][4] = 1.01974;
  phiscale[14][4] = 0.99334;
  phiscale[15][4] = 0.99997;
  phiscale[16][4] = 0.990929;
  phiscale[17][4] = 0.975634;
  phiscale[18][4] = 1.00303;
  phiscale[19][4] = 0.985434;
  phiscale[20][4] = 0.967616;
  phiscale[21][4] = 0.980996;
  phiscale[22][4] = 0.990579;
  phiscale[23][4] = 0.973296;
  phiscale[24][4] = 0.961434;
  phiscale[25][4] = 0.955421;
  phiscale[26][4] = 0.948587;
  phiscale[27][4] = 0.931421;
  phiscale[28][4] = 0.931238;
  phiscale[29][4] = 0.94423;
  phiscale[30][4] = 0.932764;
  phiscale[31][4] = 0.966975;
  phiscale[32][4] = 0.972268;
  phiscale[33][4] = 0.964495;
  phiscale[34][4] = 0.962898;
  phiscale[35][4] = 0.951517;
  phiscale[36][4] = 0.954945;
  phiscale[37][4] = 0.970193;
  phiscale[38][4] = 0.969978;
  phiscale[39][4] = 0.973643;
  phiscale[40][4] = 0.990848;
  phiscale[41][4] = 1.01816;
  phiscale[42][4] = 0.963459;
  phiscale[43][4] = 0.978044;
  phiscale[44][4] = 0.968011;
  phiscale[45][4] = 1.00723;
  phiscale[46][4] = 0.989955;
  phiscale[47][4] = 0.977922;
  phiscale[48][4] = 0.981623;
  phiscale[49][4] = 0.987167;
  phiscale[50][4] = 0.980652;
  phiscale[51][4] = 0.988063;
  phiscale[52][4] = 1.01634;
  phiscale[53][4] = 0.99809;
  phiscale[54][4] = 0.985051;
  phiscale[55][4] = 0.982958;
  phiscale[56][4] = 1.01437;
  phiscale[57][4] = 1.01895;
  phiscale[58][4] = 1.03569;
  phiscale[59][4] = 1.0199;
  phiscale[60][4] = 1.00772;
  phiscale[61][4] = 1.00965;
  phiscale[62][4] = 0.999778;
  phiscale[63][4] = 1.03219;
  phiscale[64][4] = 1.01908;
  phiscale[65][4] = 0.994986;
  phiscale[66][4] = 1.01092;
  phiscale[67][4] = 1.01101;
  phiscale[68][4] = 1.05256;
  phiscale[69][4] = 0.995279;
  phiscale[70][4] = 1.01966;
  phiscale[71][4] = 1.01925;
  phiscale[72][4] = 1.03437;
  phiscale[73][4] = 1.09619;
  phiscale[74][4] = 1.09045;
  phiscale[75][4] = 1.06572;
  phiscale[76][4] = 1.04789;
  phiscale[77][4] = 1.03944;
  phiscale[78][4] = 1.03925;
  phiscale[79][4] = 1.00293;
  phiscale[80][4] = 1.00435;
  phiscale[81][4] = 1.0026;
  phiscale[82][4] = 1.00019;
  phiscale[83][4] = 1.01719;
  phiscale[84][4] = 1.01446;
  phiscale[85][4] = 0.992784;
  phiscale[86][4] = 1.00798;
  phiscale[87][4] = 1.02701;
  phiscale[88][4] = 1.03072;
  phiscale[89][4] = 1.01656;
  phiscale[0][5] = 0.985575;
  phiscale[1][5] = 1.0037;
  phiscale[2][5] = 1.02172;
  phiscale[3][5] = 0.997933;
  phiscale[4][5] = 1.01535;
  phiscale[5][5] = 1.00884;
  phiscale[6][5] = 1.03411;
  phiscale[7][5] = 1.00496;
  phiscale[8][5] = 1.00105;
  phiscale[9][5] = 1.01106;
  phiscale[10][5] = 0.979743;
  phiscale[11][5] = 0.985568;
  phiscale[12][5] = 1.0002;
  phiscale[13][5] = 1.01733;
  phiscale[14][5] = 1.00003;
  phiscale[15][5] = 0.999219;
  phiscale[16][5] = 0.992979;
  phiscale[17][5] = 0.989186;
  phiscale[18][5] = 1.00088;
  phiscale[19][5] = 1.00454;
  phiscale[20][5] = 0.97909;
  phiscale[21][5] = 0.978088;
  phiscale[22][5] = 0.985971;
  phiscale[23][5] = 0.973762;
  phiscale[24][5] = 0.949348;
  phiscale[25][5] = 0.934134;
  phiscale[26][5] = 0.944567;
  phiscale[27][5] = 0.923076;
  phiscale[28][5] = 0.917105;
  phiscale[29][5] = 0.920272;
  phiscale[30][5] = 0.92333;
  phiscale[31][5] = 0.956102;
  phiscale[32][5] = 0.961894;
  phiscale[33][5] = 0.943972;
  phiscale[34][5] = 0.957506;
  phiscale[35][5] = 0.937969;
  phiscale[36][5] = 0.935287;
  phiscale[37][5] = 0.950902;
  phiscale[38][5] = 0.948535;
  phiscale[39][5] = 0.963971;
  phiscale[40][5] = 0.961423;
  phiscale[41][5] = 0.984611;
  phiscale[42][5] = 0.961474;
  phiscale[43][5] = 0.954815;
  phiscale[44][5] = 0.958598;
  phiscale[45][5] = 0.981048;
  phiscale[46][5] = 0.974096;
  phiscale[47][5] = 0.964914;
  phiscale[48][5] = 0.985434;
  phiscale[49][5] = 0.986682;
  phiscale[50][5] = 0.993842;
  phiscale[51][5] = 0.998059;
  phiscale[52][5] = 1.01306;
  phiscale[53][5] = 1.0039;
  phiscale[54][5] = 0.995036;
  phiscale[55][5] = 0.997697;
  phiscale[56][5] = 1.04012;
  phiscale[57][5] = 1.02815;
  phiscale[58][5] = 1.03207;
  phiscale[59][5] = 1.03419;
  phiscale[60][5] = 1.0336;
  phiscale[61][5] = 1.02819;
  phiscale[62][5] = 1.01877;
  phiscale[63][5] = 1.02579;
  phiscale[64][5] = 1.04178;
  phiscale[65][5] = 1.01904;
  phiscale[66][5] = 1.00612;
  phiscale[67][5] = 1.03288;
  phiscale[68][5] = 1.0659;
  phiscale[69][5] = 1.03144;
  phiscale[70][5] = 1.03946;
  phiscale[71][5] = 1.04802;
  phiscale[72][5] = 1.04185;
  phiscale[73][5] = 1.09881;
  phiscale[74][5] = 1.08464;
  phiscale[75][5] = 1.06695;
  phiscale[76][5] = 1.05261;
  phiscale[77][5] = 1.06486;
  phiscale[78][5] = 1.02747;
  phiscale[79][5] = 1.0059;
  phiscale[80][5] = 1.0115;
  phiscale[81][5] = 1.00673;
  phiscale[82][5] = 0.999097;
  phiscale[83][5] = 1.01395;
  phiscale[84][5] = 1.00323;
  phiscale[85][5] = 0.974085;
  phiscale[86][5] = 1.00493;
  phiscale[87][5] = 1.00697;
  phiscale[88][5] = 1.01774;
  phiscale[89][5] = 0.997998;
  phiscale[0][6] = 1.01471;
  phiscale[1][6] = 1.01996;
  phiscale[2][6] = 1.03725;
  phiscale[3][6] = 1.01556;
  phiscale[4][6] = 1.02688;
  phiscale[5][6] = 1.03699;
  phiscale[6][6] = 1.05697;
  phiscale[7][6] = 1.03617;
  phiscale[8][6] = 1.03593;
  phiscale[9][6] = 1.03186;
  phiscale[10][6] = 1.00999;
  phiscale[11][6] = 1.00247;
  phiscale[12][6] = 1.01596;
  phiscale[13][6] = 1.04507;
  phiscale[14][6] = 1.00565;
  phiscale[15][6] = 1.019;
  phiscale[16][6] = 1.02325;
  phiscale[17][6] = 1.00152;
  phiscale[18][6] = 1.00893;
  phiscale[19][6] = 1.01877;
  phiscale[20][6] = 0.992098;
  phiscale[21][6] = 1.00742;
  phiscale[22][6] = 1.01073;
  phiscale[23][6] = 0.99212;
  phiscale[24][6] = 0.961298;
  phiscale[25][6] = 0.958569;
  phiscale[26][6] = 0.945965;
  phiscale[27][6] = 0.942494;
  phiscale[28][6] = 0.930617;
  phiscale[29][6] = 0.938922;
  phiscale[30][6] = 0.94332;
  phiscale[31][6] = 0.973096;
  phiscale[32][6] = 0.972037;
  phiscale[33][6] = 0.963285;
  phiscale[34][6] = 0.960942;
  phiscale[35][6] = 0.94626;
  phiscale[36][6] = 0.941139;
  phiscale[37][6] = 0.967772;
  phiscale[38][6] = 0.966452;
  phiscale[39][6] = 0.96831;
  phiscale[40][6] = 0.990099;
  phiscale[41][6] = 1.00696;
  phiscale[42][6] = 0.967729;
  phiscale[43][6] = 0.961833;
  phiscale[44][6] = 0.967463;
  phiscale[45][6] = 0.996506;
  phiscale[46][6] = 0.994502;
  phiscale[47][6] = 0.992386;
  phiscale[48][6] = 1.00698;
  phiscale[49][6] = 0.995382;
  phiscale[50][6] = 1.00189;
  phiscale[51][6] = 1.01977;
  phiscale[52][6] = 1.03526;
  phiscale[53][6] = 1.0348;
  phiscale[54][6] = 1.02229;
  phiscale[55][6] = 1.04087;
  phiscale[56][6] = 1.06339;
  phiscale[57][6] = 1.06221;
  phiscale[58][6] = 1.06827;
  phiscale[59][6] = 1.06693;
  phiscale[60][6] = 1.07054;
  phiscale[61][6] = 1.05611;
  phiscale[62][6] = 1.06746;
  phiscale[63][6] = 1.06917;
  phiscale[64][6] = 1.072;
  phiscale[65][6] = 1.04985;
  phiscale[66][6] = 1.07572;
  phiscale[67][6] = 1.07072;
  phiscale[68][6] = 1.11823;
  phiscale[69][6] = 1.06588;
  phiscale[70][6] = 1.08315;
  phiscale[71][6] = 1.07247;
  phiscale[72][6] = 1.09258;
  phiscale[73][6] = 1.13452;
  phiscale[74][6] = 1.13953;
  phiscale[75][6] = 1.10556;
  phiscale[76][6] = 1.09693;
  phiscale[77][6] = 1.07637;
  phiscale[78][6] = 1.04735;
  phiscale[79][6] = 1.0177;
  phiscale[80][6] = 1.02167;
  phiscale[81][6] = 1.03752;
  phiscale[82][6] = 1.01901;
  phiscale[83][6] = 1.0243;
  phiscale[84][6] = 1.01699;
  phiscale[85][6] = 1.01214;
  phiscale[86][6] = 1.02332;
  phiscale[87][6] = 1.02937;
  phiscale[88][6] = 1.04012;
  phiscale[89][6] = 1.03144;
  phiscale[0][7] = 1.01467;
  phiscale[1][7] = 1.02519;
  phiscale[2][7] = 1.02364;
  phiscale[3][7] = 1.01296;
  phiscale[4][7] = 1.02752;
  phiscale[5][7] = 1.05142;
  phiscale[6][7] = 1.07096;
  phiscale[7][7] = 1.04909;
  phiscale[8][7] = 1.0444;
  phiscale[9][7] = 1.03019;
  phiscale[10][7] = 0.99943;
  phiscale[11][7] = 1.01506;
  phiscale[12][7] = 1.03175;
  phiscale[13][7] = 1.05554;
  phiscale[14][7] = 1.02356;
  phiscale[15][7] = 1.03008;
  phiscale[16][7] = 1.03975;
  phiscale[17][7] = 1.019;
  phiscale[18][7] = 1.03263;
  phiscale[19][7] = 1.02497;
  phiscale[20][7] = 1.00799;
  phiscale[21][7] = 1.01169;
  phiscale[22][7] = 1.01553;
  phiscale[23][7] = 0.994438;
  phiscale[24][7] = 0.975107;
  phiscale[25][7] = 0.967083;
  phiscale[26][7] = 0.953068;
  phiscale[27][7] = 0.936187;
  phiscale[28][7] = 0.925262;
  phiscale[29][7] = 0.935268;
  phiscale[30][7] = 0.930957;
  phiscale[31][7] = 0.954406;
  phiscale[32][7] = 0.979131;
  phiscale[33][7] = 0.962799;
  phiscale[34][7] = 0.948096;
  phiscale[35][7] = 0.952067;
  phiscale[36][7] = 0.94447;
  phiscale[37][7] = 0.968694;
  phiscale[38][7] = 0.975955;
  phiscale[39][7] = 0.974872;
  phiscale[40][7] = 0.982102;
  phiscale[41][7] = 0.996727;
  phiscale[42][7] = 0.964465;
  phiscale[43][7] = 0.965677;
  phiscale[44][7] = 0.97997;
  phiscale[45][7] = 0.97564;
  phiscale[46][7] = 1.00738;
  phiscale[47][7] = 0.986663;
  phiscale[48][7] = 0.996997;
  phiscale[49][7] = 0.998841;
  phiscale[50][7] = 1.01345;
  phiscale[51][7] = 1.03654;
  phiscale[52][7] = 1.04154;
  phiscale[53][7] = 1.04384;
  phiscale[54][7] = 1.01902;
  phiscale[55][7] = 1.03546;
  phiscale[56][7] = 1.07188;
  phiscale[57][7] = 1.0606;
  phiscale[58][7] = 1.08792;
  phiscale[59][7] = 1.08619;
  phiscale[60][7] = 1.09165;
  phiscale[61][7] = 1.08745;
  phiscale[62][7] = 1.06877;
  phiscale[63][7] = 1.08549;
  phiscale[64][7] = 1.09091;
  phiscale[65][7] = 1.07015;
  phiscale[66][7] = 1.07957;
  phiscale[67][7] = 1.09463;
  phiscale[68][7] = 1.13589;
  phiscale[69][7] = 1.09404;
  phiscale[70][7] = 1.09218;
  phiscale[71][7] = 1.10412;
  phiscale[72][7] = 1.08409;
  phiscale[73][7] = 1.13448;
  phiscale[74][7] = 1.13168;
  phiscale[75][7] = 1.12357;
  phiscale[76][7] = 1.09064;
  phiscale[77][7] = 1.0737;
  phiscale[78][7] = 1.04817;
  phiscale[79][7] = 1.03165;
  phiscale[80][7] = 1.03569;
  phiscale[81][7] = 1.03195;
  phiscale[82][7] = 1.02932;
  phiscale[83][7] = 1.03025;
  phiscale[84][7] = 1.01949;
  phiscale[85][7] = 1.00358;
  phiscale[86][7] = 1.02351;
  phiscale[87][7] = 1.02915;
  phiscale[88][7] = 1.04255;
  phiscale[89][7] = 1.02875;
  phiscale[0][8] = 0.996445;
  phiscale[1][8] = 1.01121;
  phiscale[2][8] = 1.01974;
  phiscale[3][8] = 1.0078;
  phiscale[4][8] = 1.02683;
  phiscale[5][8] = 1.02572;
  phiscale[6][8] = 1.05058;
  phiscale[7][8] = 1.03665;
  phiscale[8][8] = 1.02295;
  phiscale[9][8] = 1.02485;
  phiscale[10][8] = 0.976794;
  phiscale[11][8] = 1.00358;
  phiscale[12][8] = 1.01513;
  phiscale[13][8] = 1.05142;
  phiscale[14][8] = 1.02156;
  phiscale[15][8] = 1.02336;
  phiscale[16][8] = 1.01437;
  phiscale[17][8] = 1.00454;
  phiscale[18][8] = 1.02706;
  phiscale[19][8] = 1.013;
  phiscale[20][8] = 0.991442;
  phiscale[21][8] = 1.00127;
  phiscale[22][8] = 1.00048;
  phiscale[23][8] = 0.997361;
  phiscale[24][8] = 0.965656;
  phiscale[25][8] = 0.953285;
  phiscale[26][8] = 0.938809;
  phiscale[27][8] = 0.920781;
  phiscale[28][8] = 0.914665;
  phiscale[29][8] = 0.925177;
  phiscale[30][8] = 0.917872;
  phiscale[31][8] = 0.941332;
  phiscale[32][8] = 0.971573;
  phiscale[33][8] = 0.942704;
  phiscale[34][8] = 0.935808;
  phiscale[35][8] = 0.936559;
  phiscale[36][8] = 0.9311;
  phiscale[37][8] = 0.94884;
  phiscale[38][8] = 0.953997;
  phiscale[39][8] = 0.96696;
  phiscale[40][8] = 0.976698;
  phiscale[41][8] = 1.00033;
  phiscale[42][8] = 0.965874;
  phiscale[43][8] = 0.945433;
  phiscale[44][8] = 0.965221;
  phiscale[45][8] = 0.980194;
  phiscale[46][8] = 0.987697;
  phiscale[47][8] = 0.975757;
  phiscale[48][8] = 0.981488;
  phiscale[49][8] = 0.991588;
  phiscale[50][8] = 1.01192;
  phiscale[51][8] = 1.01238;
  phiscale[52][8] = 1.03087;
  phiscale[53][8] = 1.0358;
  phiscale[54][8] = 1.03526;
  phiscale[55][8] = 1.02556;
  phiscale[56][8] = 1.06894;
  phiscale[57][8] = 1.05098;
  phiscale[58][8] = 1.07849;
  phiscale[59][8] = 1.08654;
  phiscale[60][8] = 1.07666;
  phiscale[61][8] = 1.08323;
  phiscale[62][8] = 1.07004;
  phiscale[63][8] = 1.08665;
  phiscale[64][8] = 1.08882;
  phiscale[65][8] = 1.07795;
  phiscale[66][8] = 1.08001;
  phiscale[67][8] = 1.09727;
  phiscale[68][8] = 1.12566;
  phiscale[69][8] = 1.10629;
  phiscale[70][8] = 1.10504;
  phiscale[71][8] = 1.10664;
  phiscale[72][8] = 1.08806;
  phiscale[73][8] = 1.13582;
  phiscale[74][8] = 1.12904;
  phiscale[75][8] = 1.09942;
  phiscale[76][8] = 1.09065;
  phiscale[77][8] = 1.07385;
  phiscale[78][8] = 1.04148;
  phiscale[79][8] = 1.03787;
  phiscale[80][8] = 1.02314;
  phiscale[81][8] = 1.02636;
  phiscale[82][8] = 1.01954;
  phiscale[83][8] = 1.01695;
  phiscale[84][8] = 1.01557;
  phiscale[85][8] = 0.996516;
  phiscale[86][8] = 1.01491;
  phiscale[87][8] = 1.01614;
  phiscale[88][8] = 1.01721;
  phiscale[89][8] = 1.00725;
  phiscale[0][9] = 1.03675;
  phiscale[1][9] = 1.06331;
  phiscale[2][9] = 1.07249;
  phiscale[3][9] = 1.06956;
  phiscale[4][9] = 1.08957;
  phiscale[5][9] = 1.09853;
  phiscale[6][9] = 1.12313;
  phiscale[7][9] = 1.09689;
  phiscale[8][9] = 1.09034;
  phiscale[9][9] = 1.08967;
  phiscale[10][9] = 1.05888;
  phiscale[11][9] = 1.075;
  phiscale[12][9] = 1.08026;
  phiscale[13][9] = 1.10467;
  phiscale[14][9] = 1.08024;
  phiscale[15][9] = 1.0788;
  phiscale[16][9] = 1.09205;
  phiscale[17][9] = 1.07191;
  phiscale[18][9] = 1.08666;
  phiscale[19][9] = 1.07876;
  phiscale[20][9] = 1.06291;
  phiscale[21][9] = 1.06201;
  phiscale[22][9] = 1.06787;
  phiscale[23][9] = 1.04817;
  phiscale[24][9] = 1.02627;
  phiscale[25][9] = 1.00947;
  phiscale[26][9] = 0.983158;
  phiscale[27][9] = 0.960978;
  phiscale[28][9] = 0.956794;
  phiscale[29][9] = 0.967322;
  phiscale[30][9] = 0.958718;
  phiscale[31][9] = 0.980905;
  phiscale[32][9] = 1.00847;
  phiscale[33][9] = 0.989393;
  phiscale[34][9] = 0.97702;
  phiscale[35][9] = 0.967818;
  phiscale[36][9] = 0.967727;
  phiscale[37][9] = 0.997554;
  phiscale[38][9] = 0.996214;
  phiscale[39][9] = 0.992615;
  phiscale[40][9] = 1.01092;
  phiscale[41][9] = 1.02801;
  phiscale[42][9] = 1.01618;
  phiscale[43][9] = 0.997995;
  phiscale[44][9] = 1.00528;
  phiscale[45][9] = 1.03431;
  phiscale[46][9] = 1.03312;
  phiscale[47][9] = 1.03545;
  phiscale[48][9] = 1.03147;
  phiscale[49][9] = 1.0594;
  phiscale[50][9] = 1.06632;
  phiscale[51][9] = 1.0769;
  phiscale[52][9] = 1.10079;
  phiscale[53][9] = 1.09779;
  phiscale[54][9] = 1.07696;
  phiscale[55][9] = 1.08731;
  phiscale[56][9] = 1.12765;
  phiscale[57][9] = 1.14931;
  phiscale[58][9] = 1.15241;
  phiscale[59][9] = 1.14574;
  phiscale[60][9] = 1.15783;
  phiscale[61][9] = 1.13575;
  phiscale[62][9] = 1.1441;
  phiscale[63][9] = 1.15868;
  phiscale[64][9] = 1.15699;
  phiscale[65][9] = 1.1662;
  phiscale[66][9] = 1.16363;
  phiscale[67][9] = 1.16139;
  phiscale[68][9] = 1.16964;
  phiscale[69][9] = 1.13337;
  phiscale[70][9] = 1.1375;
  phiscale[71][9] = 1.13396;
  phiscale[72][9] = 1.13185;
  phiscale[73][9] = 1.17096;
  phiscale[74][9] = 1.16443;
  phiscale[75][9] = 1.13317;
  phiscale[76][9] = 1.11916;
  phiscale[77][9] = 1.1034;
  phiscale[78][9] = 1.07807;
  phiscale[79][9] = 1.05979;
  phiscale[80][9] = 1.05844;
  phiscale[81][9] = 1.05484;
  phiscale[82][9] = 1.04058;
  phiscale[83][9] = 1.04371;
  phiscale[84][9] = 1.03006;
  phiscale[85][9] = 1.01116;
  phiscale[86][9] = 1.03856;
  phiscale[87][9] = 1.03029;
  phiscale[88][9] = 1.05045;
  phiscale[89][9] = 1.02628;
  phiscale[0][10] = 1.04424;
  phiscale[1][10] = 1.07466;
  phiscale[2][10] = 1.05007;
  phiscale[3][10] = 1.10509;
  phiscale[4][10] = 1.10386;
  phiscale[5][10] = 1.09387;
  phiscale[6][10] = 1.13106;
  phiscale[7][10] = 1.11179;
  phiscale[8][10] = 1.11561;
  phiscale[9][10] = 1.10581;
  phiscale[10][10] = 1.07953;
  phiscale[11][10] = 1.07614;
  phiscale[12][10] = 1.08909;
  phiscale[13][10] = 1.11665;
  phiscale[14][10] = 1.0978;
  phiscale[15][10] = 1.09238;
  phiscale[16][10] = 1.10589;
  phiscale[17][10] = 1.05714;
  phiscale[18][10] = 1.11401;
  phiscale[19][10] = 1.08893;
  phiscale[20][10] = 1.07986;
  phiscale[21][10] = 1.08182;
  phiscale[22][10] = 1.09882;
  phiscale[23][10] = 1.04778;
  phiscale[24][10] = 1.03519;
  phiscale[25][10] = 1.00337;
  phiscale[26][10] = 0.99314;
  phiscale[27][10] = 0.977345;
  phiscale[28][10] = 0.961774;
  phiscale[29][10] = 0.980311;
  phiscale[30][10] = 0.959965;
  phiscale[31][10] = 1.0059;
  phiscale[32][10] = 1.00627;
  phiscale[33][10] = 1.00004;
  phiscale[34][10] = 0.988332;
  phiscale[35][10] = 0.983197;
  phiscale[36][10] = 0.969531;
  phiscale[37][10] = 0.990029;
  phiscale[38][10] = 0.989501;
  phiscale[39][10] = 1.00183;
  phiscale[40][10] = 1.00347;
  phiscale[41][10] = 1.02696;
  phiscale[42][10] = 1.00247;
  phiscale[43][10] = 1.00199;
  phiscale[44][10] = 1.00308;
  phiscale[45][10] = 1.03133;
  phiscale[46][10] = 1.05059;
  phiscale[47][10] = 1.02422;
  phiscale[48][10] = 1.0469;
  phiscale[49][10] = 1.07972;
  phiscale[50][10] = 1.055;
  phiscale[51][10] = 1.07016;
  phiscale[52][10] = 1.11415;
  phiscale[53][10] = 1.0959;
  phiscale[54][10] = 1.10447;
  phiscale[55][10] = 1.11305;
  phiscale[56][10] = 1.11762;
  phiscale[57][10] = 1.15518;
  phiscale[58][10] = 1.15056;
  phiscale[59][10] = 1.16164;
  phiscale[60][10] = 1.1635;
  phiscale[61][10] = 1.16829;
  phiscale[62][10] = 1.16344;
  phiscale[63][10] = 1.1809;
  phiscale[64][10] = 1.17163;
  phiscale[65][10] = 1.18134;
  phiscale[66][10] = 1.17977;
  phiscale[67][10] = 1.18773;
  phiscale[68][10] = 1.21138;
  phiscale[69][10] = 1.14878;
  phiscale[70][10] = 1.15335;
  phiscale[71][10] = 1.13754;
  phiscale[72][10] = 1.16469;
  phiscale[73][10] = 1.14545;
  phiscale[74][10] = 1.14571;
  phiscale[75][10] = 1.15641;
  phiscale[76][10] = 1.12408;
  phiscale[77][10] = 1.08829;
  phiscale[78][10] = 1.08828;
  phiscale[79][10] = 1.06775;
  phiscale[80][10] = 1.07567;
  phiscale[81][10] = 1.03724;
  phiscale[82][10] = 1.06441;
  phiscale[83][10] = 1.07373;
  phiscale[84][10] = 1.03672;
  phiscale[85][10] = 1.03836;
  phiscale[86][10] = 1.03597;
  phiscale[87][10] = 1.031;
  phiscale[88][10] = 1.05114;
  phiscale[89][10] = 1.02585;

  Int_t runIndex = RunIndex(run);
  Int_t phiIndex = PhiIndex(phi);

  return phiscale[phiIndex][runIndex];
}

//Momentum Correction (positive charge)
//V3.0 07/05/2004 SG 
Double_t MomCorrPos(Double_t momentum) {

  Double_t cpar[10] = {-8.49351, 41.6292, -80.363, 85.317, -55.5938, 23.1596, 
		       -6.19471, 1.03134, -0.0975565, 0.00401972};

  Double_t corr = 0.; //Additive correction
  if(momentum >= 2.5) return corr;

  double zz = cpar[0] + cpar[1]*momentum;
  for (Int_t j=2; j<10; j++) 
    zz += cpar[j]*pow(momentum, j);
  
  corr = -zz;
  return corr;
}

//Momentum Correction (negative charge)
//V3.0 07/05/2004 SG 
Double_t MomCorrNeg(Double_t momentum) {
  
  Double_t corr = 0.;
  if(momentum >= 2.0) return corr;

  Double_t zz = 7813.58*TMath::Freq(1.42733*(momentum+2.14262)) - 
    7813.67 + 0.0202517*momentum;

  corr = -zz;
  return corr;
}

//Return corrected dedx
//
// Input:

// dedx    --> dedx hit-level calibrated
// hitdedx --> number of hits used in dedx
//
// NOTE: dedx and hitdedx are usually retrieved from CdfTrack (5.X data):
//  Double_t dedx    = trk->dedxCT();
//  Int_t    hitdedx = trk->numCTHitsDedx();
//
// run     --> run number
// phi     --> track phi0 (rad)
// eta     --> track pseudorapidity
// mom     --> track momentum 
// charge  --> track charge
// 
// Output: the method return the corrected dedx
//
Double_t correctedDedx(Double_t dedx, Int_t hitdedx, Int_t run, 
		       Double_t phi, Double_t eta, Double_t mom, 
		       Int_t charge) {

  if (dedx == 0) return 0.;

  Double_t c1 = 1.;//multiplicative
  Double_t c2 = 1.;//multiplicative
  Double_t c3 = 1.;//multiplicative
  Double_t c4 = 0.;//additive

  if (charge > 0) {
    c1 = PhiCorrPos(run, phi);
    c2 = HitsCorrPos(hitdedx);
    c3 = EtaCorrPos(run, eta);
    c4 = MomCorrPos(mom);
  }
  else {
    c1 = PhiCorrNeg(run, phi);
    c2 = HitsCorrNeg(hitdedx);
    c3 = EtaCorrNeg(run, eta);
    c4 = MomCorrNeg(mom);
  }

  //Rescale COT pressure to nominal value for runs > 160457
  Double_t pcorr = 1.;
  if (run > 160457) pcorr = 1.0189286373;

  Double_t temp = dedx*c1*c2*c3*pcorr;
  temp += c4;
  return temp;
}

//Define Z = Ln(dedx/expected)
Double_t Zvar(Double_t dedx, Double_t expected) {

  if (dedx <= 0 || expected <= 0) return -999.;
  return TMath::Log(dedx) - TMath::Log(expected);
}

//Sigma z = Ln(dEdx/expected) (positive particles)
//V5.0 07/13/2004 SG (Matthew model + momentum dependence)
Double_t sigmaZPos(Int_t numused, Double_t expected, Double_t p) {

  if (numused <= 0) return -999.;

  Double_t xm = p;
  if (p >= 8.) xm = 8.0;

  Double_t scalef = 3.68166*TMath::Exp(-4.51396*xm) + 1.00204 - 0.0118334*xm;

  Double_t n0     = 61.4214;
  Double_t sigma  = 0.095909;
  Double_t b      = -1.22419;
  Double_t c      = -0.000913;
  
  Double_t n  = static_cast<Double_t>(numused); //Number of dE/dx hits
  Double_t uc = expected;                       //Predicted dE/dx
  
  return scalef*((sigma+c*(uc-15.0))*pow(n/n0,b));
}

//Sigma dEdx (positive particles)
//V5.0 07/13/2004 SG (Matthew model + momentum dependence)
Double_t sigmaDedxPos(Double_t dedx, Int_t numused, Double_t expected, Double_t p) {

  if (numused <= 0) return -999.;

  Double_t sigz = sigmaZPos(numused, expected, p);
  Double_t z    = Zvar(dedx, expected);

  return expected*TMath::Exp(z)*sigz;
}

//Sigma z = Ln(dEdx/expected) (negative particles)
//V5.0 07/13/2004 SG (Matthew model + momentum dependence)
Double_t sigmaZNeg(Int_t numused, Double_t expected, Double_t p) {

  if (numused <= 0) return -999.;

  Double_t xm = p;
  if (p >= 8.) xm = 8.0;

  Double_t scalef = 7.22507*TMath::Exp(-4.71358*xm) + 1.00122 - 0.0205669*xm;

  Double_t n0     = 60.7507;
  Double_t sigma  = 0.103933;
  Double_t b      = -0.891335;
  Double_t c      = -0.001344;
  
  Double_t n  = static_cast<Double_t>(numused); //Number of dE/dx hits
  Double_t uc = expected;                       //Predicted dE/dx
  
  return scalef*((sigma+c*(uc-15.0))*pow(n/n0,b));
}

//Sigma dEdx (negative particles)
//V5.0 07/13/2004 SG (Matthew model + momentum dependence)
Double_t sigmaDedxNeg(Double_t dedx, Int_t numused, Double_t expected, Double_t p) {

  Double_t sigz = sigmaZNeg(numused, expected, p);
  Double_t z    = Zvar(dedx, expected);

  return expected*TMath::Exp(z)*sigz;
}

//Particle dependent UC corrections
Double_t corrUC(Double_t mass, Double_t p, Int_t charge) {

    Double_t P0 = 0.;
    Double_t P1 = 0.;

    if (mass == 0.) return 0.;
    if (p == 0.) return 0.;

    if (mass < 0.053) {//Electrons
      if (charge > 0) {
	if (p <= 1.5) {
          P0 =  2.94054e-02;
          P1 = -5.83554e-02;
	}
	else if (p > 1.5 && p <= 2.5) {
          P0 = -1.80435e-01;
          P1 =  9.69496e-02;
	}
	else if (p > 2.5 && p < 8.) {
          P0 =  1.73305e-01;
          P1 = -2.84023e-02;
	}
	else {
          P0 = -0.0046;
	}
      }
      else {
        if (p <= 1.) {
          P0 =  3.16560e-01;
          P1 = -3.62277e-01;
	}
	else if (p > 2. && p < 12.) {
          P0 = -6.40744e-02;
          P1 =  1.12760e-02;
	}
	else {
          P0 = 0.071;
	}
      }
    }
    else if (mass >= 0.053 && mass < 0.122614) {//Muons
      if (charge > 0) {
	if (p <= 2.) {
          P0 = 0.019;
	}
	else if (p > 2. && p < 14.) {
          P0 =  5.35195e-02;
          P1 = -1.74744e-02;
	}
        else {
          P0 = -0.098;
	}
      }
      else {
	if (p < 2.) {
          P0 = -0.196;
	}
	else if (p >= 2. && p < 3.5) {
          P0 = -4.88814e-01;
          P1 =  1.46283e-01;
	}
	else if (p >= 3.5) {
	  P0 = -0.013;
	}
      }
    }
    else if (mass >= 0.122614 && mass < 0.317) {//Pions
      if (charge > 0) {
        if (p <= 2.) {
          P0 = -8.90625e-02;
          P1 =  1.13363e-01;
	}
        else if (p > 2. && p < 14.) {
          P0 =  1.38060e-01;
          P1 = -1.32970e-02;
	}
      }
      else {
	if (p <= 2.) {
          P0 = -4.76847e-01;
          P1 =  3.64993e-01;
	}
        else if (p > 2. && p < 14.) {
          P0 = 8.11936e-02;
          P1 = 5.78098e-03;
	}
      }
    }
    else if (mass >= 0.317 && mass < 0.716) {//Kaons
      if (charge > 0) {
        if ( p <= 2.) {
          P0 = -8.90625e-02;
          P1 =  1.13363e-01;
	}
        else if (p > 2. && p < 14.) {
          P0 = -1.57925e-01;
          P1 =  8.83067e-03;
	}
      }
      else {
        if ( p <= 2.) {
          P0 = -4.76847e-01;
          P1 =  3.64993e-01;
	}
        else if (p > 2. && p < 14.) {
          P0 = -1.96337e-01;
          P1 =  3.05811e-02;
	}
      }
    }
    else {//P corr
      if (charge > 0) {
	if (p < 0.9) {
	  P0 = 2.59e-02;
	}
	else if (p >= 0.9 && p < 4.) {
          P0 =  4.48183e-02;
          P1 = -2.10130e-02;
	}
	else {
	  P0 = -0.04;
	}
      }
      else {
	if (p < 0.7) {
	  P0 = 0.283;
	}
	else if (p >= 0.7 && p < 1.4) {
          P0 =  6.05591e-01;
          P1 = -4.61198e-01;
	}
	else if (p >= 1.4 && p < 5.) {
          P0 = -8.16198e-02;
          P1 =  3.40332e-02;
	}
	else {
	  P0 = 0.0885;
	}
      }
    }

    return P0 + P1*p;
}

//Predicted dedx (function of betagamma)
//V5.0 07/05/2004 SG (latest corrections)
//
// Input:
//  mass      --> mass hypothesis (== 0. particle specific corrections not used)
//  betagamma --> track betagamma
//  charge    --> track charge
//
Double_t predictedDedx(Double_t mass, Double_t betagamma, Int_t charge) {

  if(betagamma<=0 || charge==0) return  -999;

  Double_t x = betagamma;
  Double_t t = x/sqrt(1+x*x);  //beta

  if(x<=0) return -999;

  //Common part
  Double_t  c0   = charge>0 ?     20.1302 :    20.4296;
  Double_t  c1   = charge>0 ?      2.12314:     2.19519;
  Double_t  b    = charge>0 ?    101.26   :   105.878;
  Double_t  a1   = charge>0 ?     -4.37323:    -5.70865;
  Double_t  a2   = charge>0 ?     16.5832 :    16.3159;
  Double_t  bias = 0.;

  if (mass > 0.) {//Don't use particle dependent corrections if mass == 0.

    Double_t p = betagamma*mass;
    bias = corrUC(mass, p, charge);

  }
  
  Double_t dedx = (c1*log(x/(x+b)) + c0)/(t*t) + a1*(t-1) + a2*(t-1)*(t-1);
  return dedx + bias;
}

//Predicted dedx (function of p and mass)
//V5.0 07/05/2004 SG (latest corrections)
//
// Input:
//  mass   --> mass hypothesis (must be > 0)
//  p      --> track momentum
//  charge --> track charge
//
Double_t predictedDedxVSP(Double_t mass, Double_t p, Int_t charge) {

  if(mass <=0 || p<=0 || charge==0) return  -999;

  Double_t x = p/mass;         //betagamma

  return predictedDedx(mass, x, charge);
}

//dEdx resolution (i.e. sigma(Z) = sigma(Ln(dE/dx/expected)))
//
// Input:
//  hitdedx --> number of used hits by dedx
//  charge  --> track charge 
//  mom     --> track momentum
//  mass    --> mass hypothesis
//
// Output:
//  return sigma(Z)
//
Double_t sigmaZ(Int_t hitdedx, Int_t charge, Double_t mom, Double_t mass) {

  if (hitdedx <= 0) return -999.;

  Double_t sigma = 0.;
  if (mass == 0) {
    mass = 0.139570;
  }

  Double_t expected = predictedDedx(mass, mom, charge);

  if (charge > 0) 
    sigma = sigmaZPos(hitdedx, expected, mom);
  else
    sigma = sigmaZNeg(hitdedx, expected, mom);

  return sigma;
}

//dEdx resolution (i.e. sigma(dE/dx))
//
// Input:
//  dedx    --> measured dedx
//  hitdedx --> number of used hits by dedx
//  charge  --> track charge 
//  mom     --> track momentum
//  mass    --> mass hypothesis
//
// Output:
//  return sigma(dedx)
//
Double_t sigmaDedx(Double_t dedx, Int_t hitdedx, Int_t charge, Double_t mom, Double_t mass) {

  if (hitdedx <= 0) return -999.;

  Double_t sigma = 0.;
  if (mass == 0) {
    mass = 0.139570;
  }

  Double_t expected = predictedDedx(mass, mom, charge);

  if (charge > 0) 
    sigma = sigmaDedxPos(dedx, hitdedx, expected, mom);
  else
    sigma = sigmaDedxNeg(dedx, hitdedx, expected, mom);

  return sigma;
}

//Return corrected Z = ln(dedx/expected)
Double_t correctedZ(Double_t rawdedx, Int_t hitdedx, Int_t run, 
		    Double_t phi, Double_t eta, Double_t mom, 
		    Int_t charge, Double_t mass) {
  
  Double_t dedx     = correctedDedx(rawdedx, hitdedx, run, phi, eta, mom, charge);
  Double_t expected = predictedDedx(mass, mom, charge);
  if (dedx <= 0 || expected <= 0) return -999.;
  return Zvar(dedx,expected);
}

