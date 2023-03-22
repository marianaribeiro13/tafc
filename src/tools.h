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

  static double Norm(vector<double>);
  static vector<string> Read_File(string);
  static TSpline3* Interpolate_Photon_Spectrum(string);
  ~tools();



};

#endif
