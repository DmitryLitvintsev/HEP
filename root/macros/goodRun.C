#include <vector>
#include <algorithm>
#include <iostream>
#include "TSystem.h"

bool goodRun(const int& run) {
	
	std::vector<int> good;
	//Access the DB, make the temporary run list and save them in a STL vector
	if (run < 0) {  
		if (gSystem->AccessPathName("/cdf/atom/home/litvinse/root_macros/goodrun.sql")) {
			std::cout << "ERROR sql query file goodrun.sql not accesible !!!!" << std::endl;
			return false;
		}

		gSystem->Exec("rm -f out_tmp");
		gSystem->Exec("rm -f goodrun.list");
		gSystem->Exec("sqlplus -S cdf_reader/reader@cdfofprd @/cdf/atom/home/litvinse/root_macros/goodrun > out_tmp");
		gSystem->Exec("cat out_tmp | grep -e 1 -e 2 | awk '{print $1\" \"$2\" \"$3}'> goodrun.list");
		gSystem->Exec("rm -f out_tmp");
		
		FILE *ifp = fopen("goodrun.list","r");
		int tmp;
		Float_t xs, xe;
		Float_t online_lumi;
		while (fscanf(ifp, "%i %f %f", &tmp, &xs, &xe) != -1) {
//			online_lumi += (xe - xs);
			online_lumi += (xs);
			good.push_back(tmp);
		}
		online_lumi /= 1000.;
		std::cout << "Good run list contain: " << good.size() << " runs" << std::endl;
		std::cout << "Total online integrated luminosity: " << online_lumi 
		     << " pb^-1" << std::endl;
		fclose(ifp);
		return true;
	}
	return std::binary_search(good.begin(), good.end(), run);
}
