#ifndef Modes
#define Modes

#include "tools.h"
#include "Particle.h"
#include "Tracker.h"
#include "TApplication.h"
#include <string>
#include <chrono>
#include <fstream>
#include <string>
#include <sstream>
#include <thread>
#include "TTree.h"

void Simulation_Mode(TGeoManager* geom, int seed, int& Nmuons_total, int& Nphotons_total, 
                     int& Nphotons_detected, int& Nphotons_absorbed, int& Nphotons_lost);

void Draw_Mode(TGeoManager* geom, int seed, int N_photons_draw);

void AngularAcceptance_Mode(TGeoManager* geom, int seed, double& initial_x_muon, 
                            double& initial_y_muon, double& theta, TTree* tree);

void GeomEfficiency_Mode(TGeoManager* geom, int seed, int& Nmuons_total);

void Heatmap_Mode(TGeoManager* geom, int seed, double& phi, double& z, TTree* tree, int scintillator);

void DiskEfficiency_Mode(TGeoManager* geom, int seed, double& initial_x_muon, 
                         double& initial_y_muon, double& detector_efficiency, TTree* tree);

void DetectionEfficiency_Mode(TGeoManager* geom, int seed, int& Nphotons_total, 
                              int& Nphotons_detected, TTree* tree, int scintillator);

    

#endif