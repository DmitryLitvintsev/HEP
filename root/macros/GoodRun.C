#include "GoodRun.hh"
#include <algorithm>
#include <TSystem.h>
#include <TROOT.h>
#include <TBenchmark.h>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

GoodRun::GoodRun()
  : _initialized(true), _prevRun(-1), _prevGood(true)
{

  //Access the DB, make the temporary run list and save them in a STL vector
  if (gSystem->AccessPathName("/cdf/atom/home/litvinse/root_macros/goodrun.sql")) {
    cout << "ERROR sql query file goodrun.sql not accesible !!!!" << endl;
    _initialized=false;
  }
  
  if (_initialized) {
    gSystem->Exec("rm -f out_tmp");
    gSystem->Exec("rm -f goodrun.list");
    gSystem->Exec("sqlplus cdf_reader/reader@cdfofprd @/cdf/atom/home/litvinse/root_macros/goodrun > out_tmp");
    gSystem->Exec("cat out_tmp | grep -e 1 -e 2 | grep -v Production | grep -v Oracle | grep -v selected | awk '{print $1\" \"$2\" \"$3}'> goodrun.list");
    gSystem->Exec("rm -f out_tmp");

    FILE *ifp = fopen("goodrun.list","r");
    Int_t tmp;
    Int_t nrun_offl   = 0;
    Int_t nrun_onll   = 0;
    Int_t nrun_ononly = 0;
    Int_t nrun_wol    = 0;
    Float_t xs, xo;
    Float_t offline_lumi = 0.;
    Float_t online_lumi  = 0.;
    Float_t final_lumi   = 0.;
    while (fscanf(ifp, "%i %f %f", &tmp, &xs, &xo) != -1) {
      if (xs > 0) { 
	offline_lumi += xs;
	final_lumi   += xs;
	nrun_offl++;
      }
      if (xo > 0) { 
	online_lumi += xo;
	nrun_onll++;
      }

     if (xs <= 0 && xo > 0) { 
	final_lumi += xo;
	nrun_ononly++;
      }
      if (xs <= 0 && xo <= 0) {
	nrun_wol++;
      }
      _good.push_back(tmp);
    }
    std::sort(_good.begin(),_good.end());
    offline_lumi /= 1000.;
    online_lumi /= 1000.;
    final_lumi /= 1000.;
    cout << "Good run list contains:                   " << _good.size() << " runs" << endl;
    cout << "# runs with offline luminosity info:      " << nrun_offl << endl;
    cout << "# runs with online  luminosity info:      " << nrun_onll << endl;
    cout << "# runs with only online luminosity info:  " << nrun_ononly << endl;
    cout << "# runs without any luminosity info:       " << nrun_wol << endl;
    cout << "Offline integrated luminosity:            " << offline_lumi << endl;
    cout << "Online  integrated luminosity:            " << online_lumi << endl;
    cout << " " << endl;
    cout << "Total integrated luminosity (OFF OR ON (IF OFF 0): " << final_lumi*1.019 
         << " +/- " << final_lumi*0.06 << " pb^-1" << endl;
    cout << "NOTE: +1.9% correction applied end error set to 6% for Total lumi" << endl;
    fclose(ifp);

    //FILE *ofp = fopen("goodrun.h","w");
    ofstream ofp ("goodrun.h",ios::out);
    ofp << "\/\/Good run list contains:      " << _good.size() << " runs\n";
    ofp << "\/\/Total integrated luminosity: " << final_lumi*1.019 
        << " +/- " << final_lumi*0.06 << " pb^-1\n";
    ofp << "  Int_t goodRun[";
    ofp << _good.size() ;
    ofp << "] = {\n";
    for (Int_t irun = 0; irun <  _good.size()-1; irun++) {
       ofp << "    " << _good[irun] << ",\n";
     }
    ofp << "    " << _good[_good.size()-1] << "\n" << "  };\n";

    ofp << "  std::vector<int> goodRunVec;\n";
    ofp << "  for (Int_t rentry=0; rentry<" << _good.size() << "; rentry++) {\n";
    ofp << "      goodRunVec.push_back(goodRun[rentry]);\n";
    ofp << "  }\n";
    ofp << "  std::sort(goodRunVec.begin(),goodRunVec.end());\n";
    ofp.close();
  }
}

GoodRun::~GoodRun(){}

bool GoodRun::isGood(int run) {
  // Check cached results for speed
  if (run==_prevRun) {
    return _prevGood;
  }
  
  //Binary search on the good run STL vector
  if (_initialized) {
    _prevRun = run;
    _prevGood = std::binary_search(_good.begin(), _good.end(), run);
    return _prevGood;
  }
  else {
    cout << "There were problems making good run list.  You're SOL." << endl;
    return false;
  }
}
