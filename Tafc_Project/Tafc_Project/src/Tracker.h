#ifndef __Tracker__
#define __Tracker__

#include "Geometry.h"
#include "Particle.h"
#include "Generator.h"
#include "tools.h"
#include <vector>
#include "TGeoTrack.h"
#include "TGeoNavigator.h"
#include "TF1.h"
#include "TCanvas.h"
#include <cmath>
#include "TH2.h"
#include "Parameters.h"
using namespace std;

class Tracker
{
public:

  /////////// Constructor and Destructor //////////////////
  Tracker(TGeoManager* GeoM, Generator* gen, Particle* part);
  ~Tracker();

  //////////////Geometry check///////////
  double Update_Energy(double);
  bool CheckSameLocation();
  double FresnelLaw(double thetai, double n1, double n2);
  bool CheckReflection(double thetai, double n1, double n2);
  vector<double> GetNormal();
  //bool ReflectionHandler(int); //Removed
  bool VacuumToPlastic(double);
  bool VacuumToAluminium(double);
  bool DetectionCheck(double);
  bool Is_Photon_Detected(double E);

  /////Muon Propagators/////
  void Propagate_Muon();
  void Muon_Vacuum_Step();
  void Muon_Scintillator_Step();
  void Muon_Aluminium_Step();

  /////Photon Propagators/////
  void Propagate_Photons(int iphoton, int fphoton);
  void InitializePhotonTrack(int);
  //void Photon_Scintillator_Step(int); //Removed
  void Photon_Scintillator_Reflection_Check(int);
  // void Photon_Absorbtion(int,double); //Removed
  //void Photon_Vacuum_Step(int); //Removed
  bool Photon_Vacuum_Reflection_Check(int);
  //void Photon_Aluminium_Step(int); //Removed
  void Update_Photon(int i);

  //Data getters
  int GetN_photons(){return N_photons;};
  int GetN_photons1(){return N_photons1;};
  int GetN_absorbed(){return N_absorbed;};
  int GetN_detected(){return N_detected;};
  int GetN_lost(){return N_lost;};
  Particle* GetMuon(){return Muon;};
  bool GetDoubleCross(){return DoubleCross;};
  TGeoNavigator* GetNavigator(){return nav;};

  double CheckDensity(); //Move to Geometry
	double GetRefractiveIndex(); // Move to Geometry
	bool Check_Symmetric_Detector();
  bool Is_Detector_Region();

  double BetheBloch(double v);

  std::vector<std::pair<double,double>> PropagatePhotons_To_FirstBoundary(int iphoton, int fphoton);

  //void Reset();

private:

  TGeoManager* geom;
  TGeoNavigator* nav;
  TVirtualGeoTrack* main_track;
  double stepvalue; //defined arbitrary step
  Generator* generator; //pointer to generator object (generates random variables)
  Particle* Muon;
  vector<Particle*> Photons;
  const double *cpoint; // pointer to current position (fCurrentPoint of TGeoNavigator)
  const double *cdir; // pointer to current position (fCurrentDir of TGeoNavigator)

  int N_photons; //number of emitted photons
  int N_photons1; //number of photons 
  int N_absorbed; //number of absorbed photons in the scintillator
  int N_detected; //number of detected photons (using SIPMS)
  //int N_detected1; //number of detected photons in the first scintillator
  //int N_detected2; //number of detected photons in the second scintillator
  int N_lost; //number of photons absorbed in the aluminium
  bool DoubleCross; // This Flag Checks if the photons have been propagated or not
  bool Photons_flag; //This Flag Checks if the muon has crossed both scintillators

  Parameters param;
};

#endif