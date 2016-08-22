#include "Rtypes.h"


class CdfTrack { 
public:
	struct CdfTrackData { 
		Double_t  p;
		Double_t pt;
		Int_t charge;
		Int_t numCTHitsDedx;
		Double_t dedxCT;
		Double_t eta;
		Double_t phi0;
		Double_t tof;
		Double_t tofError;
	};
	
	CdfTrackData data;
	
	const Double_t& p()              const { return data.p; } 
	const Double_t& pt()             const { return data.pt; } 
	const Double_t& pseudoRapidity() const { return data.eta; } 
	const Double_t& phi0()           const { return data.phi0; }
	const Int_t&   charge()          const { return data.charge; }
	const Int_t&   numCTHitsDedx()   const { return data.numCTHitsDedx; }
	const Double_t& dedxCT()         const { return data.dedxCT; } 
	const Double_t& tof()            const { return data.tof; } 
	const Double_t& tofError()       const { return data.tofError; } 
	
	
	CdfTrack() {};
	~CdfTrack() {};
	
};
