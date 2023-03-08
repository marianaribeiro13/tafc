#ifndef __muon__
#define __muon__
#include <cmath>
#include "Geometry.h"

using namespace std;
class muon
{
public:
  muon(double energy, double* position, double *direction);
  ~muon();
  double GetEnergy(){return energy;};
  double GetMass(){return mass;};
  double GetAngle(){return angle;};
  double* GetMomentum(){return momentum;};
  double* GetPosition(){return position;};
  bool


protected:
  double energy;
  double angle;
  double mass = 105.6583755; //MeV
  double momentum[3];
  double position[3];



};

#endif
