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
	Geometry(double d = 0.);
	virtual ~Geometry(){std::cout<< "Destruct Geometry" << std::endl;};

protected:
    
    TGeoManager* geom; 
	TGeoVolume* top;
    double distance; //Separation between scintillators
};

#endif