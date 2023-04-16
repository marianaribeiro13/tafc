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

  /////////////// Getters ////////////////////
  double Get_stepvalue(){return stepvalue;};
  TGeoNavigator* Get_Navigator(){return nav;};
  Particle* Get_Muon(){return Muon;};
  int GetN_photons(){return N_photons;};
  int GetN_photons1(){return N_photons1;};
  int GetN_absorbed(){return N_absorbed;};
  int GetN_detected(){return N_detected;};
  int GetN_lost(){return N_lost;};
  bool GetDoubleCross(){return DoubleCross;};
  
  ///////////////// Geometry check//////////////
  bool CheckSameLocation(); //Checks whether the next position is in the same volume as the current one
  bool VacuumToPlastic(double); //Check if the muon is in vacuum near the scintillator boundary
  bool VacuumToAluminium(double); //Check if the muon is in vacuum near aluminium foil
  vector<double> GetNormal(); //Get normal vector to the next crossing boundary
  double CheckDensity(); //Get density of current material where the paricle is located
	double GetRefractiveIndex(); //Get refractive index of current material where the paricle is located
  

  //////////////Muon Propagators////////////
  void Propagate_Muon(); //Main muon propagation function - checks in which volume is the muon and calls the step functions
  void Muon_Vacuum_Step(); //Takes muon step in the vacuum (air) to the next boundary
  void Muon_Scintillator_Step(); //Propagates muon in the scintillator and calls energy loss function
  void Muon_Aluminium_Step(); //Takes muon step in the aluminium foil to the next boundary
  double Update_Energy(double); //Calls BetheBlock and updates energy and momentum of the muon 


  //////Photon Propagators//////
  void Propagate_Photons(int iphoton, int fphoton); //Propagate photons of the photon vector between the specified indexes
  void InitializePhotonTrack(int); //Initialize photon track for draw mode
  void Photon_Scintillator_Reflection_Check(int); //Computes new propagation direction after reflection or transmission in scintillator
  bool Photon_Vacuum_Reflection_Check(int); ////Computes new propagation direction after reflection or transmission in vacuum (vacuum-scintillator or vacuum-aluminium foil interfaces)
  bool CheckReflection(double thetai, double n1, double n2); //Check if photon is reflected or transmited in some boundary 
  void Update_Photon(int i); //Update photon position and track
  
  //Propagate photons until they reach the lateral disk boundary of a scintillator 
  //for the first time and save the height and phi angle positions
  std::vector<std::pair<double,double>> PropagatePhotons_To_FirstBoundary(int iphoton, int fphoton); 


  ///////Photon Detection////////////
	bool Check_Symmetric_Detector(); //Faster method for detection check (considers SIMPs at symmetric positions)
  bool Is_Detector_Region(); //Slower method for detection check (any angle for SIMps can be defined)
  bool DetectionCheck(double E); //Calls Check_Symmetric_Detector and applies SIPM detection efficiency(depends on photon energy)
  bool Is_Photon_Detected(double E); //Calls Is_Detector_Region and applies SIPM detection efficiency(depends on photon energy)


  /////////////////Applied Physics Formulas///////////////////
  double BetheBloch(double v); //Computes the average muon energy loss per unit length
  double FresnelLaw(double thetai, double n1, double n2); //Computes probability of reflection for the photon
  

private:

  double stepvalue; //defined arbitrary step

  TGeoManager* geom; //pointer to TGeoManager object (interface with Geometry)
  TGeoNavigator* nav; //pointer to TGeoNavigator object (handles propagation)
  TVirtualGeoTrack* main_track; //pointer to the main track(muon track)
  
  Generator* generator; //pointer to generator object (generates random variables)
  Particle* Muon; //pointer to muon
  vector<Particle*> Photons; //vector of photons generated during muon propagation

  const double *cpoint; // pointer to current position (fCurrentPoint of TGeoNavigator)
  const double *cdir; // pointer to current position (fCurrentDir of TGeoNavigator)

  int N_photons; //total number of propagated photons
  int N_photons1; //number of propagated photons in the top scintillator 
  int N_absorbed; //number of absorbed photons in the scintillator
  int N_detected; //number of detected photons (using SIPMS)
  int N_lost; //number of photons lost in the aluminium foil

  bool DoubleCross; //This flag checks if the muon has crossed both scintillators
  bool Photons_flag; //This flag checks if the photons have been propagated or not
  Parameters param; //struct with important parameters (defined in Parameters.h)

};

#endif
