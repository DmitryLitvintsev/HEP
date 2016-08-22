#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>
#include <cmath>
using namespace std;

int getNum(int high) {
	return (int)(( (double)rand() / (double)RAND_MAX )*(double)high);
}


template <class T>
bool readdata(const char* txt, vector<T>& runs) { 
	ifstream file(txt);
	istream_iterator<T> i_start(file),i_end;
	copy(i_start,i_end,back_inserter(runs));
	file.close();
	return true;
}

template <class T>
bool writedata(const char* txt, const vector<T>& runs) { 
	ofstream file(txt);
	copy(runs.begin(),runs.end(), ostream_iterator<T>(file,"\n"));
	file.close();
	return true;
}


