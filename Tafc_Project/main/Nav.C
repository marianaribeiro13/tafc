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
//#include <filesystem>

int main(int argc, char* argv[]){
  
    auto start = std::chrono::high_resolution_clock::now(); //Program start time
    //std::string path = std::filesystem::current_path();

    int t = time(0);

    //unsigned int nthreads = std::thread::hardware_concurrency();
    //std::cout << "Available number of threads: " << nthreads << '\n';
  
    //Initialize variables used in more than one mode
    int Nmuons_total = 0; // Total number of generated muons
    int Nphotons_total = 0; // Total number of photons (only for accepted muons)
    int Nphotons_detected = 0; // Number of photons detected in SIPMs (only for accepted muons)

    std::vector<std::thread> threads; //vector of threads (each thread has one TGeoNavigator)

    Geometry World; //Instantiate geometry object - create pointer to geometry

    Parameters param;

    //Create Muon Telescope
    World.Build_MuonTelescope(param.Radius, param.Height, param.Distance, param.Airgap, param.Thickness);

    TGeoManager* geom = World.GetGeoManager(); //Get pointer to telescope geometry
    
    geom->SetMaxThreads(param.N_threads); //Set maximum number of threads and creates thread private data for all geometry objects.
    

    ////////////////////////////////////////// SIMULATION MODE ////////////////////////////////////////////////////
    
    if(argc == 1){
      int Nphotons_absorbed = 0; // Number of photons absorbed in the scintillators
      int Nphotons_lost = 0; // Number of photons lost in the aluminium

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


    //std::string filename = path + "/ProjectResults.root";
    //TFile* outfile = new TFile(filename.c_str(),"UPDATE","ProjectResults");
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
        
        std::cout << "\n\n\n\nDiskEfficiency TTree created in ProjectResults.root\n\n";
      } else {
        tree = (TTree*)outfile->Get("DiskEfficiency");

        tree->SetBranchAddress("initial_x_muon",&initial_x_muon);
        tree->SetBranchAddress("initial_y_muon",&initial_y_muon);
        tree->SetBranchAddress("detector_efficiency",&detector_efficiency);
        outfile->Delete("DiskEfficiency;1");

        std::cout << "\n\n\n\nDiskEfficiency TTree updated in ProjectResults.root\n\n";
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

      if(!(outfile->FindKey("GeomEfficiency"))){
        tree = new TTree("GeomEfficiency","GeomEfficiency");
        
        tree->Branch("distance",&param.Distance,"distance/D");
        tree->Branch("Nmuons_total",&Nmuons_total,"Nmuons_total/I");
        tree->Branch("Nmuons_accepted",&param.Nmuons_accepted,"Nmuons_accepted/I");

        std::cout << "\n\n\n\nGeomEfficiency TTree created in ProjectResults.root\n\n";
      } else {
        tree = (TTree*)outfile->Get("GeomEfficiency");

        tree->SetBranchAddress("distance",&param.Distance);
        tree->SetBranchAddress("Nmuons_total",&Nmuons_total);
        tree->SetBranchAddress("Nmuons_accepted",&param.Nmuons_accepted);
        outfile->Delete("GeomEfficiency;1");

        std::cout << "\n\n\n\nGeomEfficiency TTree updated in ProjectResults.root\n\n";
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
    }

    ////////////////////////////////////////// DETECTOR EFFICIENCY MODE ////////////////////////////////////////////////////

    if(argc == 3 && !strcmp(argv[1],"-Det"))
    {
      int n_scintillator = 0;
      if(sscanf(argv[2],"%d",&n_scintillator)){
        std::string name = "DetectorEfficiency" + std::to_string(n_scintillator);
        
        if(!(outfile->FindKey(name.c_str()))){
          tree = new TTree(name.c_str(),name.c_str());
        
          tree->Branch("Nphotons_total",&Nphotons_total,"Nphotons_total/I");
          tree->Branch("Nphotons_detected",&Nphotons_detected,"Nphotons_detected/I");
          tree->Branch("NSIMPs",&param.n_SIPM,"NSIMPs/I");
        
        std::cout << "\n\n\n\nDetectorEfficiency" << n_scintillator << " TTree created in ProjectResults.root\n\n";
        } else {
          tree = (TTree*)outfile->Get(name.c_str());

          tree->SetBranchAddress("Nphotons_total",&Nphotons_total);
          tree->SetBranchAddress("Nphotons_detected",&Nphotons_detected);
          tree->SetBranchAddress("NSIMPs",&param.n_SIPM);
          outfile->Delete((name + ";1").c_str());

          std::cout << "\n\n\n\nDetectorEfficiency" << n_scintillator << " TTree updated in ProjectResults.root\n\n";
        }

        for (int i = 0; i < param.N_threads; i++) {
          threads.emplace_back(std::thread(DetectionEfficiency_Mode, geom, t+i, std::ref(Nphotons_total),
                                std::ref(Nphotons_detected), tree, n_scintillator));
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
        std::cout << "Number of SIPMs per scintillator: " << param.n_SIPM << '\n';
        std::cout << "Total Photons Generated: " << Nphotons_total << '\n';
        std::cout << "Photons Detected: " << Nphotons_detected << '\n';
        std::cout << "Photon Detection efficiency: " << 100*(double)Nphotons_detected/Nphotons_total << "% \n\n";
      }
    }


    ////////////////////////////////////////// HEAT MAP MODE ////////////////////////////////////////////////////

    if(argc == 2 && !strcmp(argv[1],"-Map"))
    {

      double photon_z = 0.; //Height of the photon when it hits the scintillator boundary for the first time
      double photon_phi = 0.; //Phi of the photon when it hits the scintillator boundary for the first time

      if(!(outfile->FindKey("HeatMap"))){
        tree = new TTree("HeatMap","HeatMap");
        
        tree->Branch("photon_z",&photon_z,"photon_z/D");
        tree->Branch("photon_phi",&photon_phi,"photon_phi/D");
        
        std::cout << "\n\n\n\nHeatMap TTree created in ProjectResults.root\n\n";
      } else {
        tree = (TTree*)outfile->Get("HeatMap");

        tree->SetBranchAddress("photon_z",&photon_z);
        tree->SetBranchAddress("photon_phi",&photon_phi);
        outfile->Delete("HeatMap;1");

        std::cout << "\n\n\n\nHeatMap TTree updated in ProjectResults.root\n\n";
      }

      for (int i = 0; i < param.N_threads; i++) {
        threads.emplace_back(std::thread(Heatmap_Mode, geom, t+i, std::ref(photon_phi),
                              std::ref(photon_z), tree));
      }

      //Join threads (wait for all threads to finish before continuing)
      for (auto &navigator : threads) {
        navigator.join();
      }

      outfile->Write();
      outfile->Close();
    }


    auto end = std::chrono::high_resolution_clock::now(); //Program end time
    std::chrono::duration<float> duration = end - start; //Program duration

    std::cout << "Program time: " << duration.count() << "s " << std::endl;


    return 0;
}
