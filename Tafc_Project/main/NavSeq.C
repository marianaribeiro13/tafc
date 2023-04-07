#include "Generator.h"
#include "Geometry.h"
#include "tools.h"
#include "Particle.h"
#include "Tracker.h"
#include "TApplication.h"
#include <string>
#include <chrono>


int main(int argc, char* argv[]){
  
    auto start = std::chrono::high_resolution_clock::now(); //Program start time

    double radius = 5.0;
    double height = 1.0;
    double distance = 25.;
    double airgap = .1;
    double althickness = 0.0016;
    double step = 0.001;
    int n_SIPMS = 4;
    double SIPM_size = .6;
    int N_muons = 20;
    
    Geometry World; //Instantiate geometry object - create pointer to geometry

    //Create telescope
    World.Build_MuonTelescope(radius, height, distance, airgap, althickness,n_SIPMS,SIPM_size);

    TGeoManager* geom = World.GetGeoManager(); //Get pointer to telescope geometry

    Generator* gen = new Generator();

    //Start propagation in the threads
    for (int i = 0; i < N_muons; i++) {

        Particle* Muon = gen->Generate_CosmicMuon(gen->Generate_Position(distance, height, radius));

        //Create Tracker object
        Tracker Simulation(geom, gen, Muon, step, radius, height, distance, airgap, althickness, n_SIPMS, SIPM_size);

        //Propagate muon and correspondent photons
        Simulation.Propagate_Muon();
        Simulation.Propagate_Photons(Simulation.GetN_photons());

        std::cout << '\n' << "Total Photons Generated: " << Simulation.GetN_photons() << '\n';
        std::cout << "Photons Absorbed: " << Simulation.GetN_absorbed() << '\n';
        std::cout << "Photons Detected: " << Simulation.GetN_detected() << '\n';
        std::cout << "Photons Lost: " << Simulation.GetN_lost() << '\n';

    }

    auto end = std::chrono::high_resolution_clock::now(); //Program end time
    std::chrono::duration<float> duration = end - start;
    std::cout << "Program time: " << duration.count() << "s " << std::endl;

    return 0;
}
