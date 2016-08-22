#include <ostream>
using namespace std;

class Event { 
	int run;
	int event;
public: 
	Event(const int& r, const int& e) : run(r), event(e) { } 
	~Event(){};
	Event(const Event& e) { 
		run=e.run;
		event=e.event;
	} 
	Event&  operator=(const Event& e) { 
		if ( this != &e ) { 
			run=e.run;
			event=e.event;
		}
		return *this;
	}
	bool operator==(const Event& e) const { 
		return e.event==event&&e.run==run;
	}
	bool operator!=(const Event& e) const { 
		return !operator==(e);
	}
	bool operator>(const Event& e) const { 
		return run==e.run ? event > e.event : run>e.run;
	}
	bool operator<(const Event& e) const { 
		return run==e.run ? event < e.event : run<e.run;
	}
	friend std::ostream& operator<<(std::ostream&, const Event& d);

	const int& getRun() const { return run; } 
	const int& getEvent() const { return event; } 


};

inline std::ostream&  operator<<(std::ostream& ost, const Event& d) { 
	ost << d.getRun() << " " << d.getEvent() ;
	return ost;
}
