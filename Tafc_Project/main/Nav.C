#include "Generator.h"
#include "Geometry.h"
#include "tools.h"
#include "Particle.h"
#include "Tracker.h"
#include "TApplication.h"
#include "Modes.h"
#include <string>
#include <chrono>
#include <thread>

int main(int argc, char* argv[]){
  
    auto start = std::chrono::high_resolution_clock::now(); //Program start time

    //unsigned int nthreads = std::thread::hardware_concurrency();
    //std::cout << "Available number of threads: " << nthreads << '\n';

    //////////////// PARAMETERS (all distances in centimeters) ////////////////////////
    
    double radius = 5.0; //Scintillator radius
    double height = 1.0; //Scintillator height
    double distance = 25.; //Distance between the two scintillators
    double airgap = .1; //Air gap thickness between the scintillator and the aluminium foil
    double althickness = 0.0016; //Aluminium foil thickness
    int n_SIPMS = 4; //Number of SIPMS per scintillator (not being used now)
    double SIPM_size = .6; //SIPM_size (the detector area is SIPM_size*SIPM_size)

    //Phi angle of each SIPM placed on center of lateral scintillator wall (size of the vector is the number of SIPMs per scintillator)
    std::vector<double> SIPM_angles {0., M_PI/2, M_PI, 3*M_PI/2};
    
    double step = 0.001; //Defined step for muon propagation (Energy loss in the scintillator is calculated for each step)
    
    int Nmuons_accepted = 16; // Number of muons to propagate (they have to cross both detectors)
    int N_threads = 8; // Number of threads

    ///////////////////////////////////////////////////////////////////////////////////
   
    //Initialize variables
    int Nmuons_total = 0; // Total number of generated muons
    int Nphotons_total = 0; // Total number of photons (only for accepted muons)
    int Nphotons_detected = 0; // Number of photons detected in SIPMs (only for accepted muons)
    int Nphotons_absorbed = 0; // Number of photons absorbed in the scintillators
    int Nphotons_lost = 0; // Number of photons lost in the aluminium

    std::vector<std::thread> threads; //vector of threads (each thread has one TGeoNavigator)
    int muons_per_thread = Nmuons_accepted/N_threads; // Number of accepted muons to propagate on each thread

    Geometry World; //Instantiate geometry object - create pointer to geometry

    //Create Muon Telescope
    World.Build_MuonTelescope(radius, height, distance, airgap, althickness);

    TGeoManager* geom = World.GetGeoManager(); //Get pointer to telescope geometry
    
    geom->SetMaxThreads(N_threads); //Set maximum number of threads and creates thread private data for all geometry objects.
    
    if(argc == 1){
      //Start propagation in the threads
      for (int i = 0; i < N_threads; i++) {
        threads.emplace_back(std::thread(&Modes::Simulation_Mode, geom, step, radius, height, distance, airgap, 
        althickness, n_SIPMS, SIPM_size, SIPM_angles, i, muons_per_thread,
        std::ref(Nmuons_total), std::ref(Nphotons_total), std::ref(Nphotons_detected), std::ref(Nphotons_absorbed), std::ref(Nphotons_lost)));
      }

      //Join threads (wait for all threads to finish before continuing)
      for (auto &navigator : threads) {
        navigator.join();
      }

      /////////////////////////////////////////// Print of Results /////////////////////////////////

      std::cout << "\n\n\n" << "//////////////////////// FINAL RESULTS ////////////////////////" << "\n\n";
      std::cout << "Total Muons Generated: " << Nmuons_total << '\n';
      std::cout << "Total Muons Accepted (Cross both scintillators): " << Nmuons_accepted << '\n';
      std::cout << "Geometrical acceptance: " << 100*(double)Nmuons_accepted/Nmuons_total << "% \n\n";
      std::cout << "Total Photons Generated: " << Nphotons_total << '\n';
      std::cout << "Photons Detected: " << Nphotons_detected << '\n';
      std::cout << "Photons Absorbed: " << Nphotons_absorbed << '\n';
      std::cout << "Photons Lost in Aluminium: " << Nphotons_lost << '\n';
      std::cout << "Photon Detection efficiency: " << 100*(double)Nphotons_detected/Nphotons_total << "% \n\n";
    }

    if(argc == 3 && !strcmp(argv[1],"-draw"))
    {

      int N_photons_draw = 0; //Variable that will store the number of photons drawn per muon
      if(sscanf(argv[2],"%d",&N_photons_draw)){
        //Start propagation in the threads
        for (int i = 0; i < N_threads; i++) {
          threads.emplace_back(std::thread(&Modes::Draw_Mode, geom, step, radius, height, distance, airgap, 
          althickness, n_SIPMS, SIPM_size, SIPM_angles, i, muons_per_thread, N_photons_draw));
        }

        //Join threads (wait for all threads to finish before continuing)
        for (auto &navigator : threads) {
          navigator.join();
        }

        /////////////////////////////////////////// Draw ////////////////////////////////////////
        TApplication app("app", nullptr, nullptr);
    
        TCanvas *c1 = new TCanvas("c1","c1",1200,900); //Create Canvas
    
        geom->GetTopVolume()->Draw(); //Draw geometry
        // /D: Track and first level descendents only are drawn
        // /*: Track and all descendents are drawn
        geom->DrawTracks("/*"); //Draw tracks

        app.Run();
      }
    }

    auto end = std::chrono::high_resolution_clock::now(); //Program end time
    std::chrono::duration<float> duration = end - start; //Program duration

    std::cout << "Program time: " << duration.count() << "s " << std::endl;


    return 0;
}
