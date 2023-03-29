#ifndef __Geometry__
#define __Geometry__

#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include <iostream>

class Geometry{

public:

	// constructors and destructor
	Geometry();
	virtual ~Geometry() = default;
	double GetDistance(){return Distance;};
	double GetHeight(){return Height;};
	double GetRadius(){return Radius;};

	// Build Telescope methods
	void Build_MuonTelescope(double radius, double height, double distance,  double airgap, double althickness); // radius and height of scintillator and distance between scintillators

protected:
    double Radius;
		double innerradius;
		double outerradius;
		double Height;
		double Distance;
		double MaxHeight;
		double MaxRadius;
		double Thickness;
		double Airgap;
    TGeoManager* geom; //Manager pointer to geometry
};

#endif
