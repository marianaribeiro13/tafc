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

class Modes{
public:
    Modes() = default;
    //void Heatmap_Mode(Tracker*,DataManager*);
    //void HeatmapSingle_Mode(Tracker*,DataManager*);
    //void GeomEfficiency_Mode(Tracker*,DataManager*);
    //void EMap_Mode(Tracker*,DataManager*);
    //void Draw_Mode(Tracker*, int Nphotons);
    //void Simulation_Mode(Tracker*);
    static void Draw_Mode(TGeoManager* geom, double step, double radius, double height, double distance, 
                            double airgap, double althickness, int n_SIPMS, double SIPM_size, std::vector<double> SIPM_angles, 
                            int seed, int N_muons, int N_photons_draw);

    static void Simulation_Mode(TGeoManager* geom, double step, double radius, double height, double distance, 
                            double airgap, double althickness, int n_SIPMS, double SIPM_size, std::vector<double> SIPM_angles, int seed, int N_muons,
                            int& Nmuons_total, int& Nphotons_total, int& Nphotons_detected, int& Nphotons_absorbed, int& Nphotons_lost);
    ~Modes() = default;
};

#endif