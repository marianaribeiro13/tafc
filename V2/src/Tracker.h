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
  Tracker(double,double,double,double,double,double,Generator*,int,double);
  ~Tracker();

  double Update_Energy(double);
  bool CheckSameLocation();
  double FresnelLaw(double,double,double);
  bool CheckReflection(double,double,double);
  vector<double> GetNormal();
  bool VacuumToPlastic(double);
  bool VacuumToAluminium(double);
  bool DetectionCheck(const double*,double);

  /////Simulation Control/////
  void Reset();
  void Generate_New_Muon();
  void Insert_New_Muon(Particle*);

  /////Propagators/////
  void Propagate_Muon();
  void Muon_Vacuum_Step();
  void Muon_Scintillator_Step();
  void Muon_Aluminium_Step();


  void Propagate_Photons();
  void InitializePhotonTrack(int);
  void Photon_Scintillator_Reflection_Check(int);
  bool Photon_Vacuum_Reflection_Check(int);
  void Update_Photon_Track(int);

  //Data
  int GetN_photons(){return N_photons;};
  int GetN_absorbed(){return N_absorbed;};
  int GetN_detected(){return N_detected;};
  int GetN_lost(){return N_lost;};
  Particle* GetMuon(){return Muon;};
  bool GetDoubleCross(){return DoubleCross;};
  bool GetPhotonFlag(){return Photons_flag;};


  //Draw Mode
  void Draw(int n);
  void Fill_Heatmap(TH2D*);


private:
  double stepvalue;
  int MainTrackID;
  Generator* generator;
  Particle* Muon;
  vector<Particle*> Photons;
  TF1* BetheBloch;
  const double *cpoint;
  TVirtualGeoTrack* main_track;

  int N_photons;
  int N_absorbed;
  int N_detected;
  int N_lost;

  bool DoubleCross;
  bool Photons_flag;


};

#endif
