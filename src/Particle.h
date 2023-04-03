#ifndef __Particle__
#define __Particle__
#include <cmath>
#include <vector>
#include "TObject.h"

using namespace std;

class Particle : public TObject {

public:

  //constructors and destructor
  Particle(int const PDG, double p, std::vector<double> const& d,vector<double> const& x); // give pdg, momentum and direction vector
  Particle(Particle* part); //copy constructor
  //Particle(int const PDG, double E, std::vector<double> const& d); // give pdg, energy and direction vector
  ~Particle();

  //getters
  int const GetPDG() {return pdg;};
  double const GetMass() {return mass;};
  double GetMomentum() {return momentum;};
  double GetEnergy() {return energy;};
  double GetVelocity(){return momentum/energy;};
  vector<double> GetDirection(){return direction;};
  vector<double> GetStartingPosition(){return StartingPosition;};
  double GetTime(){return time;};
  vector<double> GetPosition(){return Position;};
  //change variables
  void ChangeEnergy(double E){ energy = E; };
  void ChangeMomentum(double p){ momentum = p; };
  void ChangeDirection(std::vector<double> const& d){ direction = d; };
  void ChangePosition(const double*);

  void IncreaseTime(double dt){time += dt;};

  //Calculate variables
  double CalculateMomentum(double E){return sqrt(E*E-mass*mass);};
  double CalculateEnergy(double p){return sqrt(p*p+mass*mass);};

private:

  int pdg; //Particle PDG
  double mass; //mass in MeV
  double energy; //Energy in MeV
  double momentum; //Momentum in MeV
  std::vector<double> direction; //Direction
  vector<double> StartingPosition;
  vector<double> Position;
  double time;

};

#endif
