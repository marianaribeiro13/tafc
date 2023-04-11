#include "Generator.h"
#include "Geometry.h"
#include "tools.h"
#include "Particle.h"
#include "Tracker.h"
#include "TApplication.h"
#include "Modes.h"
#include "TTree.h"
#include "TFile.h"
#include <string>
#include <chrono>
#include <thread>

int main(int argc, char* argv[]){
  
    auto start = std::chrono::high_resolution_clock::now(); //Program start time
    
    int t = time(0);

    //unsigned int nthreads = std::thread::hardware_concurrency();
    //std::cout << "Available number of threads: " << nthreads << '\n';
  
    //Initialize variables
    int Nmuons_total = 0; // Total number of generated muons
    int Nphotons_total = 0; // Total number of photons (only for accepted muons)
    int Nphotons_detected = 0; // Number of photons detected in SIPMs (only for accepted muons)
    int Nphotons_absorbed = 0; // Number of photons absorbed in the scintillators
    int Nphotons_lost = 0; // Number of photons lost in the aluminium

    std::vector<std::thread> threads; //vector of threads (each thread has one TGeoNavigator)

    Geometry World; //Instantiate geometry object - create pointer to geometry

    Parameters param;

    //Create Muon Telescope
    World.Build_MuonTelescope(param.Radius, param.Height, param.Distance, param.Airgap, param.Thickness);

    TGeoManager* geom = World.GetGeoManager(); //Get pointer to telescope geometry
    
    geom->SetMaxThreads(param.N_threads); //Set maximum number of threads and creates thread private data for all geometry objects.
    

    ////////////////////////////////////////// SIMULATION MODE ////////////////////////////////////////////////////
    
    if(argc == 1){
      //Start propagation in the threads
      for (int i = 0; i < param.N_threads; i++) {
        threads.emplace_back(std::thread(Simulation_Mode, geom, t+i, std::ref(Nmuons_total),
                                          std::ref(Nphotons_total), std::ref(Nphotons_detected), 
                                          std::ref(Nphotons_absorbed), std::ref(Nphotons_lost)));
      }

      //Join threads (wait for all threads to finish before continuing)
      for (auto &navigator : threads) {
        navigator.join();
      }

      /////////////////////////////////////////// Print of Results /////////////////////////////////

      std::cout << "\n\n\n" << "//////////////////////// FINAL RESULTS ////////////////////////" << "\n\n";
      std::cout << "Total Muons Generated: " << Nmuons_total << '\n';
      std::cout << "Total Muons Accepted (Cross both scintillators): " << param.Nmuons_accepted << '\n';
      std::cout << "Geometrical acceptance: " << 100*(double)param.Nmuons_accepted/Nmuons_total << "% \n\n";
      std::cout << "Total Photons Generated: " << Nphotons_total << '\n';
      std::cout << "Photons Detected: " << Nphotons_detected << '\n';
      std::cout << "Photons Absorbed: " << Nphotons_absorbed << '\n';
      std::cout << "Photons Lost in Aluminium: " << Nphotons_lost << '\n';
      std::cout << "Photon Detection efficiency: " << 100*(double)Nphotons_detected/Nphotons_total << "% \n\n";
    }


    ////////////////////////////////////////// DRAW MODE ////////////////////////////////////////////////////
    
    if(argc == 3 && !strcmp(argv[1],"-draw"))
    {

      int N_photons_draw = 0; //Variable that will store the number of photons drawn per muon
      if(sscanf(argv[2],"%d",&N_photons_draw)){
        //Start propagation in the threads
        for (int i = 0; i < param.N_threads; i++) {
          threads.emplace_back(std::thread(Draw_Mode, geom, t+i, N_photons_draw));
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


    TFile* outfile = new TFile("~/Tafc_Project/ProjectResults.root","UPDATE","ProjectResults");
    TTree *tree;

    ////////////////////////////////////////// DISK EFFICIENCY MODE ////////////////////////////////////////////////////

    if(argc == 2 && !strcmp(argv[1],"-Disk"))
    {

      double initial_x_muon = 0.;
      double initial_y_muon = 0.;
      double detector_efficiency = 0.;

      if(!(outfile->FindKey("DiskEfficiency"))){
        tree = new TTree("DiskEfficiency","DiskEfficiency");
        
        tree->Branch("initial_x_muon",&initial_x_muon,"initial_x_muon/D");
        tree->Branch("initial_y_muon",&initial_y_muon,"initial_y_muon/D");
        tree->Branch("detector_efficiency",&detector_efficiency,"detector_efficiency/D");
        
        std::cout << "DiskEfficiency TTree created \n\n";
      } else {
        tree = (TTree*)outfile->Get("DiskEfficiency");

        tree->SetBranchAddress("initial_x_muon",&initial_x_muon);
        tree->SetBranchAddress("initial_y_muon",&initial_y_muon);
        tree->SetBranchAddress("detector_efficiency",&detector_efficiency);
        outfile->Delete("DiskEfficiency;1");

        std::cout << "DiskEfficiency TTree updated \n\n";
      }

      for (int i = 0; i < param.N_threads; i++) {
        threads.emplace_back(std::thread(DiskEfficiency_Mode, geom, t+i, std::ref(initial_x_muon),
                              std::ref(initial_y_muon), std::ref(detector_efficiency), tree));
      }

      //Join threads (wait for all threads to finish before continuing)
      for (auto &navigator : threads) {
        navigator.join();
      }

      outfile->Write();
      outfile->Close();

    }


    ////////////////////////////////////////// GEOMETRY ACCEPTANCE EFFICIENCY MODE ////////////////////////////////////////////////////

    if(argc == 2 && !strcmp(argv[1],"-Geo"))
    {
      //TFile* outfile = new TFile("~/Tafc_Project/GeomEfficiency.root","UPDATE","GeomEfficiency");
      //TTree *tree;

      if(!(outfile->FindKey("GeomEfficiency"))){
        tree = new TTree("GeomEfficiency","GeomEfficiency");
        
        tree->Branch("distance",&param.Distance,"distance/D");
        tree->Branch("Nmuons_total",&Nmuons_total,"Nmuons_total/I");
        tree->Branch("Nmuons_accepted",&param.Nmuons_accepted,"Nmuons_accepted/I");

        std::cout << "GeomEfficiency TTree created \n\n";
      } else {
        tree = (TTree*)outfile->Get("GeomEfficiency");

        tree->SetBranchAddress("distance",&param.Distance);
        tree->SetBranchAddress("Nmuons_total",&Nmuons_total);
        tree->SetBranchAddress("Nmuons_accepted",&param.Nmuons_accepted);
        outfile->Delete("GeomEfficiency;1");

        std::cout << "GeomEfficiency TTree updated \n\n";
      }

      for (int i = 0; i < param.N_threads; i++){
        threads.emplace_back(std::thread(GeomEfficiency_Mode, geom, t+i, std::ref(Nmuons_total)));
      }

      //Join threads (wait for all threads to finish before continuing)
      for (auto &navigator : threads) {
        navigator.join();
      }

      tree->Fill();
      outfile->Write();
      outfile->Close();

      /////////////////////////////////////////// Print of Results /////////////////////////////////

      std::cout << "\n\n\n" << "//////////////////////// FINAL RESULTS ////////////////////////" << "\n\n";
      std::cout << "Distance between scintillators: " << param.Distance << '\n';
      std::cout << "Total Muons Generated: " << Nmuons_total << '\n';
      std::cout << "Total Muons Accepted (Cross both scintillators): " << param.Nmuons_accepted << '\n';
      std::cout << "Geometrical acceptance: " << 100*(double)param.Nmuons_accepted/Nmuons_total << "% \n\n";

      std::cout << "GeomEfficiency TTree added to ProjectResults.root \n\n";
    }


    auto end = std::chrono::high_resolution_clock::now(); //Program end time
    std::chrono::duration<float> duration = end - start; //Program duration

    std::cout << "Program time: " << duration.count() << "s " << std::endl;


    return 0;
}
