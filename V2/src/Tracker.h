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
using namespace std;

class Tracker : public Geometry
{
public:
  Tracker(double,double,double,double,double,double,Generator*);
  ~Tracker();
  double CheckDensity();
  double Update_Energy(double);
  bool CheckSameLocation();
  double FresnelLaw(double,double,double);
  bool CheckReflection(double,double,double);
  vector<double> GetNormal();
  bool ReflectionHandler(int);
  bool VacuumToPlastic(double);
  bool VacuumToAluminium(double);

  /////Propagators/////
  void Propagate_Muon();
  void Muon_Vacuum_Step();
  void Muon_Scintillator_Step();
  void Muon_Aluminium_Step();


  void Propagate_Photons();
  void InitializePhotonTrack(int);
  void Photon_Vacuum_Step(int);
  void Photon_Scintillator_Step(int);
  void Photon_Scintillator_Absorbtion(int,double);
  void Photon_Aluminium_Step(int);

  //Data
  int GetN_photons(){return N_photons;};
  int GetN_absorbed(){return N_absorbed;};
  int GetN_detected(){return N_detected;};
  int GetN_lost(){return N_lost;};

  //Draw Mode
  void Draw();
  void Propagate_Photons_DrawMode(int);

  //Debug Mode//
  void print_vector(const double*);
  void Debug();

private:
  double stepvalue;
  Generator* generator;
  Particle* Muon;
  vector<Particle*> Photons;
  TF1* BetheBloch;
  const double *cpoint;


  int N_photons;
  int N_absorbed;
  int N_detected;
  int N_lost;
  bool DoubleCross;


};

#endif
