#include "Generator.h"
#include "Geometry.h"
#include "tools.h"
#include "Particle.h"
#include "Tracker.h"
#include "TApplication.h"
#include <string>
#include <chrono>
#include <thread>
#include <mutex>

mutex mu;
//function for multithreading
void navigation_thread(TGeoManager* geom, double step, double radius, double height, 
double distance, double airgap, double althickness, int n_SIPMS, double SIPM_size, int i, int N_muons)
{
  Generator* gen = new Generator(i); // I think this should be inside lock but it is working fine

  for(int j; j < N_muons; j++){
    //Locked stuff happens one thread at a time
    mu.lock(); //When we create a particle we create a TDatabasePDG instance (the same in every thread - seg fault) - that is why we lock
    //Generator* gen = new Generator(i);
    //Generate random muon in the scintillator incident plane
    Particle* Muon = gen->Generate_CosmicMuon(gen->Generate_Position(distance, height, radius));
    mu.unlock();

    //Create Tracker object
    Tracker Simulation(geom, gen, Muon, step, radius, height, distance, airgap, althickness, n_SIPMS, SIPM_size);

    //Propagate muon and correspondent photons
    Simulation.Propagate_Muon();
    Simulation.Propagate_Photons(Simulation.GetN_photons());

    //Prints for current thread
    mu.lock(); 
    //std::cout << "numero navigators na thread: " << gGeoManager->GetListOfNavigators()->GetEntries() << '\n';
    std::cout << "\n\n" << "ThreadID: " << Simulation.GetNavigator()->GetThreadId();
    std::cout << '\n' << "Total Photons Generated: " << Simulation.GetN_photons() << '\n';
    std::cout << "Photons Absorbed: " << Simulation.GetN_absorbed() << '\n';
    std::cout << "Photons Detected: " << Simulation.GetN_detected() << '\n';
    std::cout << "Photons Lost: " << Simulation.GetN_lost() << '\n';
    mu.unlock();

  }

  delete gen;
}

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
    int N_navigators = 10;
    int N_muons = 20;
    int muons_per_thread = N_muons/N_navigators;

    std::vector<std::thread> navigators; //vector of threads (each thread has one TGeoNavigator)

    Geometry World; //Instantiate geometry object - create pointer to geometry

    //Create telescope
    World.Build_MuonTelescope(radius, height, distance, airgap, althickness,n_SIPMS,SIPM_size);

    TGeoManager* geom = World.GetGeoManager(); //Get pointer to telescope geometry
    
    geom->SetMaxThreads(N_navigators); //Set maximum number of threads and creates thread private data for all geometry objects.

    //Start propagation in the threads
    for (int i = 0; i < N_navigators; i++) {
      navigators.emplace_back(std::thread(navigation_thread, geom, step, radius, height, distance, airgap, althickness, n_SIPMS, SIPM_size, i, muons_per_thread));
    }

    //Join threads (wait for all threads to finish before continuing)
    for (auto &navigator : navigators) {
      navigator.join();
    }

    auto end = std::chrono::high_resolution_clock::now(); //Program end time
    std::chrono::duration<float> duration = end - start;
    std::cout << "Program time: " << duration.count() << "s " << std::endl;

    return 0;
}
