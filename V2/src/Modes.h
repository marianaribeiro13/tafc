#ifndef __Modes__
#define __Modes__

#include "tools.h"
#include "Particle.h"
#include "Tracker.h"
#include "DataManager.h"
#include "TApplication.h"
#include <string>
#include <chrono>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

void Heatmap_Mode(Tracker*,DataManager*);
void HeatmapSingle_Mode(Tracker*,DataManager*);
void EMap_Mode(Tracker*,DataManager*);
void Draw_Mode(Tracker*,DataManager*,int);
void Simulation_Mode(Tracker*,DataManager*,char*);

#endif
