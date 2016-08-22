#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>
#include <iterator>
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


class Duplicate {
public:
	Duplicate(int r=0,int e=0,double m=0, double em=0,double ptP=0, double ptp=0) : _run(r),_event(e),_mass(m), _emass(em), _ptP(ptP),_ptp(ptp) { }
	~Duplicate() { } 


	bool operator>(const Duplicate& d) const { 
		return  _run==d._run ? ( _event==d._event ? _mass>d._mass : _event>d._event ) : _run>d._run;
	}

	bool operator<(const Duplicate& d) const { 
		return  _run==d._run ? ( _event==d._event ? _mass<d._mass : _event<d._event ) : _run<d._run;
	}

	bool operator==(const Duplicate& d) const { 
		return (_run==d._run&&
			_event==d._event&&
			fabs(_mass-d._mass)<1.e-16);
	}

	bool operator!=(const Duplicate& d) const { 
		return !operator==(d); 
	}

	bool operator==(const pair<int,int>& p) const { 
		pair<int,int> tmp;
		return (tmp==p);
	}

	bool operator!=(const pair<int,int>& p) const { 
		return  !operator==(p);
	}
		
	const int&     event()   const { return _event;}
	const int&     run()     const { return _run;  }
	const double&  mass()    const { return _mass; }
	const double&  emass()   const { return _emass; }
	const double&  Pt()      const { return _ptP; }
	const double&  pt()      const { return _ptp; }

	friend std::istream& operator>>(std::istream&, Duplicate&);
	friend std::ostream& operator<<(std::ostream&, const Duplicate& d);

private:
	int _run;
	int _event;
	double _mass;
	double _emass;
	double _ptP;
	double _ptp;
};

std::istream&  operator>>(std::istream& ist, Duplicate& d) { 
	int count; 
	char txt[256];
	sprintf(txt,"%d %d %.10g %.10g %.10g %.10g",d._run,d._event,d._mass,d._emass,d._ptP,d._ptp);
		 
//	ist >> d._run >> d._event >> d._mass >> d._emass >> d._ptP >> d._ptp;
	ist >> txt;
	return ist;
}

std::ostream&  operator<<(std::ostream& ost, const Duplicate& d) { 
	ost << d.run() << " " << d.event() << " " << d.mass() <<  " " << d.emass() << " " << d.Pt() << " " << d.pt() ;
	return ost;
}

