#ifndef __tools__
#define __tools__
#include <iostream>
#include <cstdio>
#include <vector>
#include <cmath>
#include "TF1.h"
#include <fstream>
#include <string>
#include <sstream>
#include "TSpline.h"
#include "TRandom.h"

using namespace std;

class tools{
public:

  static double Norm(vector<double>);   //returns the norm of a vector
  static double Angle_Between_Vectors(vector<double>, vector<double>); //returns angle between two vectors
  static double SnellLaw(double thetai, double n1, double n2);
  static vector<double> Get_Reflected_Dir(vector<double> i, vector<double> n); //Get reflected direction vector, from incident and normal vectors
  static vector<double> Get_Refracted_Dir(vector<double> i, vector<double> n, double thetai, double n1, double n2); //Get refracted direction vector, from incident and normal directions
  static vector<string> Read_File(string); //Reads a file to a vector of strings
  static TSpline3* Interpolate_Photon_Spectrum(string); //Receives a filename and outputs an interpolation
  static vector<double> NormalizeVector(vector<double>);
  static void print_vector(const double*);
  ~tools();



};

#endif
