#ifndef __generators__
#define __generators__

#include "tools.h"
#include "Photon.h"
using namespace std;

class Generator{
public:
  Generator();
  vector<double> Generate_Vector(); //creates a random, normalized 3d vector
  double Random_Distribution(double ,double ,TF1* ); //generates a random number according to a distribution
  vector<double> Generate_Position(double); //Genertes a muon's starting position
  vector<double> Random_Distribution_2D(TF1*,double,double,double,double,double); //Generates 2 random numbers according to a pdf
  double Generate_Photon_Energy(); //Generates photon energy according to an interpolation of a pdf, use tools::Interpolate_Photon_Spectrum to get the distribution
  int Generate_Photon_Number(double); //Generates an integer according to a Poisson distribution arround the input
  //Muon* Generate_Muon(double);
  //Photon* Generate_Photon(vector<double>);
  double Uniform(double xmin,double xmax){return Random->Uniform(xmin,xmax);};
  ~Generator();

private:
  TRandom *Random;
  TF1 *Momentum_Distribution;
  TSpline3 *Photon_Spectrum;
};

#endif
