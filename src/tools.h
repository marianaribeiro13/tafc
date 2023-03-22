#ifndef __tools__
#define __tools__
#include <iostream>
#include <cstdio>
#include <vector>
#include <cmath>
#include "TF1.h"
#include "Muon.h"
#include <fstream>
#include <string>
#include <sstream>
#include "TSpline.h"
#include "TRandom.h"

using namespace std;

class tools{
public:

  static double Norm(vector<double>);   //returns the norm of a vector
  static vector<string> Read_File(string); //Reads a file to a vector of strings
  static TSpline3* Interpolate_Photon_Spectrum(string); //Receives a filename and outputs an interpolation
  ~tools();



};

#endif
