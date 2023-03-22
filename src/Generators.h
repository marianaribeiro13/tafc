#ifndef __generators__
#define __generators__

#include "tools.h"


class Generator{
public:
  Generator();
  vector<double> Generate_Vector(); //creates a random, normalized 3d vector
  double Random_Distribution(double ,double ,TF1* ); //generates a random number according to a distribution
  vector<double> Generate_Position(double,double,double);
  vector<double> Random_Distribution_2D(TF1*,double,double,double,double,double);
  double Generate_Photon_Energy(TSpline3*);
  int Generate_Photon_Number(double);

  ~Generator();

private:
  TRandom *Random;
};

#endif
