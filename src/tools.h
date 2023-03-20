#ifndef __tools__
#define __tools__
#include <iostream>
#include <cstdio>
#include <vector>
#include <cmath>
#include "TF1.h"
#include "Muon.h"

using namespace std;

class tools{
public:
  static vector<double> Generate_Vector(); //creates a random, normalized 3d vector
  static double Random_Distribution(double ,double ,TF1* ); //generates a random number according to a distribution
  static  muon* Generate_Muon(vector<double>);
  static vector<double> Generate_Position(double,double,double);
  static double Norm(vector<double>);
  static vector<double> Random_Distribution_2D(TF1*,double,double,double,double,double); /
  ~tools();

};

#endif
