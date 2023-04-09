#ifndef __Geometry__
#define __Geometry__

#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include <iostream>
#include "TGeoTrack.h"
#include "TGeoNavigator.h"
#include "tools.h"
#include <vector>
using namespace std;
class Geometry{

public:

	// constructors and destructor
	Geometry(double radius, double height, double distance,  double airgap, double althickness,int n_SIPMS,double SIPM_size);
	~Geometry();
	double GetDistance(){return Distance;};
	double GetHeight(){return Height;};
	double GetRadius(){return Radius;};
	TGeoManager* GetGeoManager(){return geom;};
	bool Check_Symmetric_Detector(const double*);
	vector<double> GetNormal(double,double,vector<double>,const double*);
	bool VacuumToPlastic(double,double);
	bool VacuumToAluminium(double,double);



	// Build Telescope methods
	void Build_MuonTelescope(double radius, double height, double distance,  double airgap, double althickness,int n_SIPMS,double SIPM_size); // radius and height of scintillator and distance between scintillators

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
	double SIPM_alpha;

  TGeoManager* geom; //Manager pointer to geometry

};

#endif
