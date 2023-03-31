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
	TGeoManager* GetGeoManager(){return geom;};
	bool CheckDetector(const double*);

	// Build Telescope methods
	void Build_MuonTelescope(double radius, double height, double distance,  double airgap, double althickness,int n_SIPMS); // radius and height of scintillator and distance between scintillators

protected:
    double Radius;
		double Height;
		double Distance;
		double innerradius;
		double outerradius;
		double Thickness;
		double Airgap;
		double MaxHeight;
		double MaxRadius;
		int n_SIPM;
		double SIPM_size;
		double SIPM_angle;
		double SIPM_angular_acceptance;

    TGeoManager* geom; //Manager pointer to geometry
};

#endif
