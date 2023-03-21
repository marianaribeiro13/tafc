#ifndef __Tracking__
#define __Tracking__

#include "Geometry.h"
#include "TVirtualGeoTrack.h"
#include <vector>
#include "Muon.h"
#include "tools.h"

class Tracking : public Geometry{

public:
	// constructors and destructor
	Tracking();
	~Tracking();

	void Propagate();

	void DefinedStep(double stepvalue);

	double CrossNextBoundary();

	void Draw();

	/*// getters
	int Ndim() const {return ndim;}
	double T() const {return x[0];}
	double X(int i) const {return x[i+1];}
	double* GetArray() {return x;}

	// operators
	ODEpoint operator*(double) const;
	ODEpoint operator+(const ODEpoint&) const;
	ODEpoint operator-(const ODEpoint&) const;

	void operator=(const ODEpoint&);

	const double& operator[] (int) const;
	double& operator[] (int);

	// print
	friend ostream& operator<<(ostream&, const ODEpoint&);*/

protected:

	double muon_step;
	TVirtualGeoTrack* track;
	Muon* muon;
};

#endif
