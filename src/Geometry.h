#ifndef __Geometry__
#define __Geometry__

#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"


class Geometry{

public:
	// constructors and destructor
	Geometry(double d = 0.);
	~Geometry();

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
    
    TGeoManager* geom; 
	TGeoVolume* top;
    double distance; //Separation between scintillators
};

#endif