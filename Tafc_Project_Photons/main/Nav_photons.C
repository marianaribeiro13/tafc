#include "Generator.h"
#include "Geometry.h"
#include "tools.h"
#include "Particle.h"
#include "Tracker.h"
#include "TApplication.h"
#include <string>
#include <chrono>
#include <thread>

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
    int N_navigators_photons = 100;

    std::vector<std::thread> navigators_photons; //vector of threads (each thread has one TGeoNavigator)

    Geometry World; //Instantiate geometry object - create pointer to geometry

    //Create telescope
    World.Build_MuonTelescope(radius, height, distance, airgap, althickness,n_SIPMS,SIPM_size);

    TGeoManager* geom = World.GetGeoManager(); //Get pointer to telescope geometry

    Generator* gen = new Generator();
    Particle* Muon = gen->Generate_CosmicMuon(gen->Generate_Position(distance, height, radius));

    Tracker* Simulation= new Tracker(geom, gen, Muon, step, radius, height, distance, airgap, althickness, n_SIPMS, SIPM_size);
    Simulation->Propagate_Muon();
    std::cout << '\n' << "Total Photons Generated: " << Simulation->GetN_photons() << '\n';
    std::cout << "Photons Absorbed: " << Simulation->GetN_absorbed() << '\n';
    std::cout << "Photons Detected: " << Simulation->GetN_detected() << '\n';
    std::cout << "Photons Lost: " << Simulation->GetN_lost() << '\n';
    Simulation->Propagate_Photons(Simulation->GetN_photons());

    auto end = std::chrono::high_resolution_clock::now(); //Program end time
    std::chrono::duration<float> duration = end - start;
    std::cout << "Program time: " << duration.count() << "s " << std::endl;

    return 0;
}
