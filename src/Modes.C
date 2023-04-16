#include "Modes.h"
#include <mutex>
#include <numeric>

std::mutex mu;


//////////////////////////////////////// SIMULATION MODE //////////////////////////////////////

void Simulation_Mode(TGeoManager* geom, int seed,
int& Nmuons_total, int& Nphotons_total, int& Nphotons_detected, int& Nphotons_absorbed, int& Nphotons_lost)
{
  std::thread::id this_id = std::this_thread::get_id(); //Thread id

  //Initialize variables
  int Nphotons_total_thisthread = 0;
  int Nphotons_detected_thisthread = 0;
  int Nphotons_absorbed_thisthread = 0;
  int Nphotons_lost_thisthread = 0;

  Generator* gen = new Generator(seed); // I think this should be inside lock but it is working fine

  Parameters param;

  int N_muons = param.muons_per_thread;

  for(int j = 0; j < N_muons; j++){

    //Generate random muon in the scintillator incident plane to add to the simulation
    Particle* Muon = gen->Generate_CosmicMuon(gen->Generate_Position(param.Distance, param.Height, param.Radius));

    //Create Tracker object
    Tracker Simulation(geom, gen, Muon);

    //Propagate muon and create emited photons
    Simulation.Propagate_Muon();

    //Check if the muon crossed both scintillators (If it did propagate photons)
    if(Simulation.GetDoubleCross()){
      Simulation.Propagate_Photons(0,Simulation.GetN_photons());
    
      Nphotons_total_thisthread += Simulation.GetN_photons();
      Nphotons_detected_thisthread += Simulation.GetN_detected();
      Nphotons_absorbed_thisthread += Simulation.GetN_absorbed();
      Nphotons_lost_thisthread += Simulation.GetN_lost();
      
    } else {
      N_muons++; //The muon was not accepted - propagate one more
    }
  }
  
  delete gen;

  //Locked stuff happens one thread at a time
  mu.lock();
  std::cout << "thread " << this_id << "\n\n"; //Print thread id
  //std::cout << "\n\n" << "ThreadID: " << Simulation.GetNavigator()->GetThreadId();
  Nmuons_total += N_muons;
  Nphotons_total += Nphotons_total_thisthread;
  Nphotons_detected += Nphotons_detected_thisthread;
  Nphotons_absorbed += Nphotons_absorbed_thisthread;
  Nphotons_lost += Nphotons_lost_thisthread;
  mu.unlock();
}


//////////////////////////////////////// DRAW MODE ///////////////////////////////////////////////

void Draw_Mode(TGeoManager* geom, int seed, int N_photons_draw)
{
  std::thread::id this_id = std::this_thread::get_id(); //Thread id

  Generator* gen = new Generator(seed); // I think this should be inside lock but it is working fine

  Parameters param;
  int N_muons = param.muons_per_thread;

  for(int j = 0; j < N_muons; j++){

    //Generate random muon in the scintillator incident plane to add to the simulation
    Particle* Muon = gen->Generate_CosmicMuon(gen->Generate_Position(param.Distance, param.Height, param.Radius));

    //Create Tracker object
    Tracker Simulation(geom, gen, Muon);

    //Propagate muon and create emited photons
    Simulation.Propagate_Muon();

    //Check if the muon crossed both scintillators (If it did propagate photons)
    if(Simulation.GetDoubleCross()){
      int m = Simulation.GetN_photons()/N_photons_draw;
     
      for(int k=0; k<N_photons_draw; k++){
        int i=m*k;
        Simulation.Propagate_Photons(i,i+1);
      }
    } else {
      N_muons++; //The muon was not accepted - propagate one more
    }
  }

  delete gen;
}


///////////////////////////////////// ANGULAR ACCEPTANCE MODE ////////////////////////////////////////

void AngularAcceptance_Mode(TGeoManager* geom, int seed, double& initial_x_muon, 
                          double& initial_y_muon, double& theta, TTree* tree){

  std::thread::id this_id = std::this_thread::get_id(); //Thread id

  Generator* gen = new Generator(seed); // I think this should be inside lock but it is working fine

  Parameters param;
  int N_muons = param.muons_per_thread;

  std::vector<double> n {0,0,-1};

  for(int j = 0; j < N_muons; j++){

    //Generate random muon in the scintillator incident plane to add to the simulation
    Particle* Muon = gen->Generate_CosmicMuon(gen->Generate_Position(param.Distance, param.Height, param.Radius));

    //Create Tracker object
    Tracker Simulation(geom, gen, Muon);

    //Propagate muon and create emited photons
    Simulation.Propagate_Muon();

    //Check if the muon crossed both scintillators (If it did propagate photons)
    if(Simulation.GetDoubleCross()){
      
      std::vector<double> d = Muon->GetDirection();
      
      mu.lock();
      initial_x_muon = Muon->GetStartingPosition()[0];
      initial_y_muon = Muon->GetStartingPosition()[1];
      theta = acos(std::inner_product(d.begin(), d.end(), n.begin(), 0.));
      tree->Fill();
      mu.unlock();

    } else {
      N_muons++; //The muon was not accepted - propagate one more
    }

  }
  
  delete gen;

  //Locked stuff happens one thread at a time
  mu.lock();
  std::cout << "thread " << this_id << "\n\n"; //Print thread id
  mu.unlock();
}


///////////////////////////////////// GEOMETRICAL EFFICIENCY MODE ///////////////////////////////////////

void GeomEfficiency_Mode(TGeoManager* geom, int seed, int& Nmuons_total)
{
  std::thread::id this_id = std::this_thread::get_id(); //Thread id

  Generator* gen = new Generator(seed); // I think this should be inside lock but it is working fine

  Parameters param;
  int N_muons = param.muons_per_thread;

  for(int j = 0; j < N_muons; j++){

    //Generate random muon in the scintillator incident plane to add to the simulation
    Particle* Muon = gen->Generate_CosmicMuon(gen->Generate_Position(param.Distance, param.Height, param.Radius));

    //Create Tracker object
    Tracker Simulation(geom, gen, Muon);

    //Propagate muon and create emited photons
    Simulation.Propagate_Muon();

    //Check if the muon crossed both scintillators (If it did propagate photons)
    if(!(Simulation.GetDoubleCross())){
      N_muons++; //The muon was not accepted - propagate one more
    }
  }
  
  delete gen;

  //Locked stuff happens one thread at a time
  mu.lock();
  std::cout << "thread " << this_id << "\n\n"; //Print thread id
  //std::cout << "\n\n" << "ThreadID: " << Simulation.GetNavigator()->GetThreadId();
  Nmuons_total += N_muons;
  mu.unlock();
}


/////////////////////////////////////////// HEATMAP MODE ////////////////////////////////////////

void Heatmap_Mode(TGeoManager* geom, int seed, double& phi, double& z, TTree* tree, int scintillator){

  std::thread::id this_id = std::this_thread::get_id(); //Thread id

  Generator* gen = new Generator(seed); // I think this should be inside lock but it is working fine

  Parameters param;
  int N_muons = param.muons_per_thread;

   std::vector<std::pair<double,double>> points;

  for(int j = 0; j < N_muons; j++){

    //Generate random muon in the scintillator incident plane to add to the simulation
    Particle* Muon = gen->Generate_CosmicMuon(gen->Generate_Position(param.Distance, param.Height, param.Radius));

    //Create Tracker object
    Tracker Simulation(geom, gen, Muon);

    //Propagate muon and create emited photons
    Simulation.Propagate_Muon();

    //Check if the muon crossed both scintillators (If it did propagate photons)
    if(Simulation.GetDoubleCross()){
      if(scintillator == 1){

        points = Simulation.PropagatePhotons_To_FirstBoundary(0,Simulation.GetN_photons1());
      
      }else{ //scintillator == 2
        
        points = Simulation.PropagatePhotons_To_FirstBoundary(Simulation.GetN_photons1(),Simulation.GetN_photons());
        
      } 
      
      
      mu.lock();
      for(auto point : points){
        phi = point.first;
        z = point.second;
        tree->Fill();
      }
      mu.unlock();

    } else {
      N_muons++; //The muon was not accepted - propagate one more
    }

  }
  
  delete gen;

  //Locked stuff happens one thread at a time
  mu.lock();
  std::cout << "thread " << this_id << "\n\n"; //Print thread id
  mu.unlock();
  
}




/////////////////////////////////////////////////// DISK EFFICIENCY MODE /////////////////////////////////////

void DiskEfficiency_Mode(TGeoManager* geom, int seed, double& initial_x_muon, 
                          double& initial_y_muon, double& detector_efficiency, TTree* tree){

  std::thread::id this_id = std::this_thread::get_id(); //Thread id

  Generator* gen = new Generator(seed); // I think this should be inside lock but it is working fine

  Parameters param;
  int N_muons = param.muons_per_thread;

  for(int j = 0; j < N_muons; j++){

    //Generate random muon in the scintillator incident plane to add to the simulation
    Particle* Muon = gen->Generate_CosmicMuon(gen->Generate_Position(param.Distance, param.Height, param.Radius));

    //Create Tracker object
    Tracker Simulation(geom, gen, Muon);

    //Propagate muon and create emited photons
    Simulation.Propagate_Muon();

    //Check if the muon crossed both scintillators (If it did propagate photons)
    if(Simulation.GetDoubleCross()){
      Simulation.Propagate_Photons(0,Simulation.GetN_photons1());
      
      mu.lock();
      detector_efficiency = 100*(double)Simulation.GetN_detected()/Simulation.GetN_photons1();
      initial_x_muon = Muon->GetStartingPosition()[0];
      initial_y_muon = Muon->GetStartingPosition()[1];
      tree->Fill();
      mu.unlock();

    } else {
      N_muons++; //The muon was not accepted - propagate one more
    }

  }
  
  delete gen;

  //Locked stuff happens one thread at a time
  mu.lock();
  std::cout << "thread " << this_id << "\n\n"; //Print thread id
  mu.unlock();
}



void DetectionEfficiency_Mode(TGeoManager* geom, int seed, int& Nphotons_total, int& Nphotons_detected, TTree* tree, int scintillator){

  std::thread::id this_id = std::this_thread::get_id(); //Thread id

  //Initialize variables
  int Nphotons_total_thisthread = 0;
  int Nphotons_detected_thisthread = 0;

  Generator* gen = new Generator(seed); // I think this should be inside lock but it is working fine

  Parameters param;

  int N_muons = param.muons_per_thread;

  for(int j = 0; j < N_muons; j++){

    //Generate random muon in the scintillator incident plane to add to the simulation
    Particle* Muon = gen->Generate_CosmicMuon(gen->Generate_Position(param.Distance, param.Height, param.Radius));

    //Create Tracker object
    Tracker Simulation(geom, gen, Muon);

    //Propagate muon and create emited photons
    Simulation.Propagate_Muon();

    //Check if the muon crossed both scintillators (If it did propagate photons)
    if(Simulation.GetDoubleCross()){

      if(scintillator == 1){

        Simulation.Propagate_Photons(0,Simulation.GetN_photons1());
        Nphotons_total_thisthread += Simulation.GetN_photons1();

      }else{ //scintillator == 2

        Simulation.Propagate_Photons(Simulation.GetN_photons1(),Simulation.GetN_photons());
        Nphotons_total_thisthread += Simulation.GetN_photons() - Simulation.GetN_photons1();
      } 
    
      Nphotons_detected_thisthread += Simulation.GetN_detected();

    } else {
      N_muons++; //The muon was not accepted - propagate one more
    }
  }
  
  delete gen;

  //Locked stuff happens one thread at a time
  mu.lock();
  std::cout << "thread " << this_id << "\n\n"; //Print thread id
  //std::cout << "\n\n" << "ThreadID: " << Simulation.GetNavigator()->GetThreadId();
  Nphotons_total += Nphotons_total_thisthread;
  Nphotons_detected += Nphotons_detected_thisthread;
  mu.unlock();
}






