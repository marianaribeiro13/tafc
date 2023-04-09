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
double distance, double airgap, double althickness, int n_SIPMS, double SIPM_size, int i, int N_muons,
int& Nmuons_total, int& Nphotons_total, int& Nphotons_detected)
{
  int Nphotons_total_thisthread = 0;
  int Nphotons_detected_thisthread = 0;

  Generator* gen = new Generator(i); // I think this should be inside lock but it is working fine

  for(int j; j < N_muons; j++){
   
    //Locked stuff happens one thread at a time
    /*mu.lock(); //When we create a particle we create a TDatabasePDG instance (the same in every thread - seg fault) - that is why we lock
    //Generate random muon in the scintillator incident plane and add it to the simulation
    Particle* Muon = gen->Generate_CosmicMuon(gen->Generate_Position(distance, height, radius));
    mu.unlock();*/

    Particle* Muon = gen->Generate_CosmicMuon(gen->Generate_Position(distance, height, radius));

    //Create Tracker object
    Tracker Simulation(geom, gen, Muon, step, radius, height, distance, airgap, althickness, n_SIPMS, SIPM_size);

    //Propagate muon and correspondent photons
    Simulation.Propagate_Muon();
    //Check if the muon crossed both scintillators (If)
    if(Simulation.GetDoubleCross()){
      Simulation.Propagate_Photons(Simulation.GetN_photons());
      Nphotons_total_thisthread += Simulation.GetN_photons();
      Nphotons_detected_thisthread += Simulation.GetN_detected();
      
    } else {
      N_muons++; //The muon was not accepted - propagate one more
    }
  }
  
  delete gen;

  //Locked stuff happens one thread at a time
  mu.lock();
  //std::cout << "\n\n" << "ThreadID: " << Simulation.GetNavigator()->GetThreadId();
  Nmuons_total += N_muons;
  Nphotons_total += Nphotons_total_thisthread;
  Nphotons_detected += Nphotons_detected_thisthread;
  mu.unlock();
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
    int N_threads = 16; // Number of threads
    int Nmuons_accepted = 1000; // Number of muons to propagate (they have to cross both detectors)
    int muons_per_thread = Nmuons_accepted/N_threads; // Number of muons to propagate on each thread

    int Nmuons_total = 0;
    int Nphotons_total = 0;
    int Nphotons_detected = 0; // photons detected in SIPMs

    std::vector<std::thread> navigators; //vector of threads (each thread has one TGeoNavigator)

    Geometry World; //Instantiate geometry object - create pointer to geometry

    //Create telescope
    World.Build_MuonTelescope(radius, height, distance, airgap, althickness,n_SIPMS,SIPM_size);

    TGeoManager* geom = World.GetGeoManager(); //Get pointer to telescope geometry
    
    geom->SetMaxThreads(N_threads); //Set maximum number of threads and creates thread private data for all geometry objects.

    //Start propagation in the threads
    for (int i = 0; i < N_threads; i++) {
      navigators.emplace_back(std::thread(navigation_thread, geom, step, radius, height, distance, airgap, 
      althickness, n_SIPMS, SIPM_size, i, muons_per_thread,
      std::ref(Nmuons_total), std::ref(Nphotons_total), std::ref(Nphotons_detected)));
    }

    //Join threads (wait for all threads to finish before continuing)
    for (auto &navigator : navigators) {
      navigator.join();
    }
 
    std::cout << "\n\n\n" << "//////////////////////// FINAL RESULTS ////////////////////////" << "\n\n";
    std::cout << "Total Muons Generated: " << Nmuons_total << '\n';
    std::cout << "Total Muons Accepted (Cross both scintillators): " << Nmuons_accepted << '\n';
    std::cout << "Geometrical acceptance: " << 100*Nmuons_accepted/Nmuons_total << "% \n\n";
    std::cout << "Total Photons Generated: " << Nphotons_total << '\n';
    std::cout << "Photons Detected: " << Nphotons_detected << '\n';
    std::cout << "Photon Detection efficiency: " << 100*Nphotons_detected/Nphotons_total << "% \n\n";


    //std::cout << "Photons Absorbed: " << Simulation.GetN_absorbed() << '\n';
    //std::cout << "Photons Lost: " << Simulation.GetN_lost() << '\n';

    auto end = std::chrono::high_resolution_clock::now(); //Program end time
    std::chrono::duration<float> duration = end - start;
    std::cout << "Program time: " << duration.count() << "s " << std::endl;

    return 0;
}
