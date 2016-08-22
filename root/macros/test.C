#include <stdio.h>
#include <stdlib.h>
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"

using namespace std;

int test(const char* name) {
	system("rm -f root.txt");
	system("ls *root > root.txt");
 	char fname[256];
 	FILE* fp = fopen("root.txt","r");
 	while (fscanf(fp,"%s",&fname)!=EOF) { 
		 if (fname[0] == '\0' || fname[0] == '#') continue;
		 TFile   f(fname);
		 TTree*  t = (TTree*)f.Get(name);
		 char txt[100];
		 sprintf(txt,"mv %s bad",fname);
		 int rc=0;
		 if (t) { 
			 if (t->GetEntries()>1.23456e+09) { 
				 rc = system(txt);
				 if (rc) { 
					 printf("failed to move file %s\n",fname);
				 }
			 }
		 }
		 else { 
				 rc = system(txt);
				 if (rc) { 
					 printf("failed to move file %s\n",fname);
				 }
		 }
	}
	return 0;
}
