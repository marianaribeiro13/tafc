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
#include "TH2.h"

class DataManager
{
public:
  DataManager();
  ~DataManager();
  void Fill_Efficiency_Map(Tracker*);
  void Draw_Efficiency_Map(string);
  void Fill_Heatmap(Tracker*);
  void Draw_Heatmap(string);
private:

  TH2D *Efficiency_Map;
  TH2D *Heatmap;

};
#endif
