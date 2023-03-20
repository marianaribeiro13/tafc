#ifndef __Muon__
#define __Muon__
#include <cmath>
#include <vector>
#include "TObject.h"

class Muon : public TObject{
public:
  Muon(std::vector<double> x, double theta, double E);
  ~Muon();
  int GetPDG() {return pdg;};
  double GetMomentum() {return momentum;};
  double GetVelocity() {return velocity;};
  double GetEnergy() {return energy;};
  std::vector<double> GetPosition();
  std::vector<double> GetDirection();
  void ChangeEnergy(const double &E){ energy = E; };
  void ChangeMomentum(){ momentum = sqrt(energy*energy - mass*mass); };
  void ChangeVelocity(){ velocity = double(momentum)/energy; };
  void ChangePosition(const std::vector<double> &x){ position = x; };
  void ChangeDirection(const std::vector<double> &d){ direction = d; };

private:
  int pdg;
  double mass = 105.6583755; //MeV
  double energy;
  double velocity;
  double momentum;
  std::vector<double> position; //Four position
  std::vector<double> direction; //Direction

};

#endif
