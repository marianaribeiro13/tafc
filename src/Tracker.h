#ifndef __Tracker__
#define __Tracker__

#include "Geometry.h"
#include "Particle.h"
#include "Generators.h"
#include "tools.h"
#include <vector>
#include "TGeoTrack.h"
#include "TGeoNavigator.h"
#include "TF1.h"
#include "TCanvas.h"
#include <cmath>
#include "TH2.h"
using namespace std;

class Tracker : public Geometry
{
public:

  /////////// Constructor and Destructor //////////////////
  Tracker(double radius, double height, double distance, double airgap, 
  double althickness, double step, Generator* g, int n_SIPMS,double SIPM_size);
  ~Tracker();

  //////////////Geometry check///////////
  double Update_Energy(double);
  bool CheckSameLocation();
  double FresnelLaw(double,double,double);
  bool CheckReflection(double,double,double);
  vector<double> GetNormal();
  //bool ReflectionHandler(int); //Removed
  bool VacuumToPlastic(double);
  bool VacuumToAluminium(double);
  bool DetectionCheck(const double*,double);

  /////Muon Propagators/////
  void Propagate_Muon();
  void Muon_Vacuum_Step();
  void Muon_Scintillator_Step();
  void Muon_Aluminium_Step();

  /////Photon Propagators/////
  void Propagate_Photons(int n);
  void InitializePhotonTrack(int);
  //void Photon_Scintillator_Step(int); //Removed
  void Photon_Scintillator_Reflection_Check(int);
  // void Photon_Absorbtion(int,double); //Removed
  //void Photon_Vacuum_Step(int); //Removed
  bool Photon_Vacuum_Reflection_Check(int);
  //void Photon_Aluminium_Step(int); //Removed
  void Update_Photon_Track(int i);

  //Data getters
  int GetN_photons(){return N_photons;};
  int GetN_absorbed(){return N_absorbed;};
  int GetN_detected(){return N_detected;};
  int GetN_lost(){return N_lost;};
  Particle* GetMuon(){return Muon;};
  bool GetDoubleCross(){return DoubleCross;};

  //Draw Mode
  void Draw();
  void Heatmap();

private:

  double stepvalue; //defined arbitrary step
  Generator* generator; //pointer to generator object (generates random variables)
  Particle* Muon;
  TVirtualGeoTrack* main_track; //Muon track (primary track)
  vector<Particle*> Photons;
  TF1* BetheBloch;
  const double *cpoint;


  int N_photons; //number of emitted photons
  int N_absorbed; //number of absorbed photons
  int N_detected; //number of detected photons (using SIPMS)
  int N_lost;
  bool DoubleCross;


};

#endif
