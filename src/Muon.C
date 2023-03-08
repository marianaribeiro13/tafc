#include "Muon.h"

muon::muon(double e, double* p, double *d)
{
  energy = e;
  angle = acos((-d[2])/sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]));
  for(int i=0;i<3;i++)
  {
    position[i]= p[i];
    momentum[i]= sqrt(e*e-mass*mass) * d[i];
  };

}

muon::~muon(){}
