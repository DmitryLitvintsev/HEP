//
// A Root macro to:
//  - retrieve from DB good runs list
//  - select good runs
//
// Usage:
//     goodRun(-1)        --> Access the DB and make the selected god run list
//                            NOTE: goodRun(-1) retrun false in case of errors
//                            NOTE2: should be called once in the job
//
//     goodRun(runnumber) --> return true for good runs and false for bad runs
//  
// V0.1 Stefano Giagu
//
#include <vector>
#include <algorithm>
#include <TSystem.h>
#include <TROOT.h>
#include <TBenchmark.h>
#include <iostream>

using namespace std;

vector<Int_t> good;

Bool_t goodRun1(Int_t run, double& lum=0) {

  //Access the DB, make the temporary run list and save them in a STL vector
  if (run < 0) {  
    if (gSystem->AccessPathName("/cdf/atom/home/litvinse/root_macros/goodrun.sql")) {
      cout << "ERROR sql query file goodrun.sql not accesible !!!!" << endl;
      return false;
    }

    gSystem->Exec("rm -f out_tmp");
    gSystem->Exec("rm -f goodrun.list");
    gSystem->Exec("sqlplus cdf_reader/reader@cdfofprd @/cdf/atom/home/litvinse/root_macros/goodrun > out_tmp");
    gSystem->Exec("cat out_tmp | grep -e 1 -e 2 | grep -v Production | grep -v Oracle | grep -v selected | awk '{print $1\" \"$2\" \"$3}'> goodrun.list");
    gSystem->Exec("rm -f out_tmp");
    
    FILE *ifp = fopen("goodrun.list","r");
    Int_t tmp;
    Float_t xs, xe;
    Float_t online_lumi;
    while (fscanf(ifp, "%i %f %f", &tmp, &xs, &xe) != -1) {
      online_lumi += (xe - xs);
      good.push_back(tmp);
    }
    online_lumi /= 1000.;
    cout << "Good run list contain: " << good.size() << " runs" << endl;
    cout << "Total online integrated luminosity: " << online_lumi 
         << " pb^-1" << endl;
    fclose(ifp);
    return true;
  }
  
  //Binary search on the good run STL vector
  return binary_search(good.begin(), good.end(), run);

}
