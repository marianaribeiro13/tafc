#ifndef __DataManager__
#define __DataManager__

#include "Tracker.h"
#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include <string>
#include <fstream>
#include <sstream>
#include "TGraph2D.h"

class DataManager
{
public:
  DataManager(int);
  ~DataManager();
  void Extract_Data(Tracker*,int);
  void Print_Data(int);
  void Draw_Efficiency_Graph();
private:
  int size;
  vector<double> x;
  vector<double> y;
  vector<double> ZenithAngle;
  vector<int> N_photons;
  vector<int> N_detected;
  vector<double> Efficiency;
  vector<int> N_absorbed;
  vector<int> N_lost;
  vector<int> DoubleCross;

};
#endif
