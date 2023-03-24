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

	// Build Telescope methods
	void Build_MuonTelescope(double radius, double height, double distance,  double airgap, double althickness); // radius and height of scintillator and distance between scintillators

protected:
    
    TGeoManager* geom; //Manager pointer to geometry
};

#endif