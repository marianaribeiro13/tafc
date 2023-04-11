#ifndef __Modes__
#define __Modes__

#include "tools.h"
#include "Particle.h"
#include "Tracker.h"
//#include "DataManager.h"
#include "TApplication.h"
#include <string>
#include <chrono>
#include <fstream>
#include <string>
#include <sstream>
#include <thread>
#include "TTree.h"

    //void Heatmap_Mode(Tracker*,DataManager*);
    //void HeatmapSingle_Mode(Tracker*,DataManager*);
    //void GeomEfficiency_Mode(Tracker*,DataManager*);
    //void EMap_Mode(Tracker*,DataManager*);
    //void Draw_Mode(Tracker*, int Nphotons);
    //void Simulation_Mode(Tracker*);

void Simulation_Mode(TGeoManager* geom, int seed, int& Nmuons_total, int& Nphotons_total, 
                     int& Nphotons_detected, int& Nphotons_absorbed, int& Nphotons_lost);

void Draw_Mode(TGeoManager* geom, int seed, int N_photons_draw);

void DiskEfficiency_Mode(TGeoManager* geom, int seed, double& initial_x_muon, 
                         double& initial_y_muon, double& detector_efficiency, TTree* tree);

void GeomEfficiency_Mode(TGeoManager* geom, int seed, int& Nmuons_total);
    

#endif