#include "Modes.h"
#include <mutex>

std::mutex mu;

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

  int N_muons = param.Nmuons_accepted;

  for(int j = 0; j < N_muons; j++){

    //Generate random muon in the scintillator incident plane to add to the simulation
    Particle* Muon = gen->Generate_CosmicMuon(gen->Generate_Position(param.Distance, param.Height, param.Radius));

    //Create Tracker object
    Tracker Simulation(geom, gen, Muon);

    //Propagate muon and create emited photons
    Simulation.Propagate_Muon();

    //Check if the muon crossed both scintillators (If it did propagate photons)
    if(Simulation.GetDoubleCross()){
      Simulation.Propagate_Photons(Simulation.GetN_photons());
    
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



void Draw_Mode(TGeoManager* geom, int seed, int N_photons_draw)
{
  std::thread::id this_id = std::this_thread::get_id(); //Thread id

  Generator* gen = new Generator(seed); // I think this should be inside lock but it is working fine

  Parameters param;
  int N_muons = param.Nmuons_accepted;

  for(int j = 0; j < N_muons; j++){

    //Generate random muon in the scintillator incident plane to add to the simulation
    Particle* Muon = gen->Generate_CosmicMuon(gen->Generate_Position(param.Distance, param.Height, param.Radius));

    //Create Tracker object
    Tracker Simulation(geom, gen, Muon);

    //Propagate muon and create emited photons
    Simulation.Propagate_Muon();

    //Check if the muon crossed both scintillators (If it did propagate photons)
    if(Simulation.GetDoubleCross()){
      Simulation.Propagate_Photons(N_photons_draw);
    } else {
      N_muons++; //The muon was not accepted - propagate one more
    }
  }
  
  delete gen;
}

void DiskEfficiency_Mode(TGeoManager* geom, int seed, double& initial_x_muon, 
                          double& initial_y_muon, double& detector_efficiency, TTree* tree){

  std::thread::id this_id = std::this_thread::get_id(); //Thread id

  Generator* gen = new Generator(seed); // I think this should be inside lock but it is working fine

  Parameters param;
  int N_muons = param.Nmuons_accepted;

  for(int j = 0; j < N_muons; j++){

    //Generate random muon in the scintillator incident plane to add to the simulation
    Particle* Muon = gen->Generate_CosmicMuon(gen->Generate_Position(param.Distance, param.Height, param.Radius));
    initial_x_muon = Muon->GetStartingPosition()[0];
    initial_y_muon = Muon->GetStartingPosition()[1];

    //Create Tracker object
    Tracker Simulation(geom, gen, Muon);

    //Propagate muon and create emited photons
    Simulation.Propagate_Muon();

    //Check if the muon crossed both scintillators (If it did propagate photons)
    if(Simulation.GetDoubleCross()){
      Simulation.Propagate_Photons(Simulation.GetN_photons());
      detector_efficiency = Simulation.GetN_detected()/Simulation.GetN_photons();
      tree->Fill();

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


void GeomEfficiency_Mode(TGeoManager* geom, int seed, int& Nmuons_total)
{
  std::thread::id this_id = std::this_thread::get_id(); //Thread id

  Generator* gen = new Generator(seed); // I think this should be inside lock but it is working fine

  Parameters param;
  int N_muons = param.Nmuons_accepted;

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



/*void Heatmap_Mode(Tracker* T, DataManager* Data)
{
  for(int i=0;i<1000;i++)
  {
    T->Propagate_Muon();
    Data->Fill_Heatmap(T);
    T->Generate_New_Muon();
  }

  Data->Draw_Heatmap("Photon_Heatmap.pdf");

  return;
}

void HeatmapSingle_Mode(Tracker* T,DataManager* Data)
{
  T->Propagate_Muon();
  Data->Fill_Heatmap(T);
  T->Generate_New_Muon();
  Data->Draw_Heatmap("Single_Photon_Heatmap.pdf");
}

void GeomEfficiency_Mode(Tracker* T,DataManager* Data)
{
  int N = 1000;
  int n=0;
  for(int i =0;i<N;i++)
  {
    cout<<i<<endl;
    T->Propagate_Muon();
    if(T->GetDoubleCross())
    {
      n++;
    }
    T->Generate_New_Muon();
  }
  cout<<(double) n/N<<endl;
  return;
}*/

// void EMap_Mode(Tracker* T, DataManager* Data)
// {
//   vector<double> d = {0,0,-1};
//   vector<double> x = {0,0,0};
//   x[2] = T->GetHeight()+T->GetDistance()/2 - 1e-12;
//
//   Particle* M0 = new Particle(13,1000,d,x);
//   T->Insert_New_Muon(M0);
//   T->Propagate_Muon();
//   T->Propagate_Photons();
//   Data->Fill_Efficiency_Map(T);
//
//   for(int i=0;i<10;i++)
//   {
//     cout<<i<<endl;
//     for(int j=0;j<10;j++)
//     {
//       x[0] = (double) (i/10) * T->GetRadius() * cos((double) (j/100) * 2 *M_PI);
//       x[1] = (double) (i/10) * T->GetRadius() * sin((double) (j/100) * 2 *M_PI);
//
//       Particle* M = new Particle(13,1000,d,x);
//       T->Insert_New_Muon(M);
//       T->Propagate_Muon();
//       T->Propagate_Photons();
//       Data->Fill_Efficiency_Map(T);
//
//     }
//   }
//
//   Data->Draw_Efficiency_Map("Efficiency_Map.pdf");
//   return;
// }

/*void Draw_Mode(Tracker T, int Nphotons)
{


  T.Propagate_Muon();
  TApplication app("app", nullptr, nullptr);
  T->Draw(n);
  app.Run();

  return;
}

void Simulation_Mode(Tracker T, int& Nphotons_total, int& Nphotons_detected, int& Nphotons_absorbed, int& Nphotons_lost)
{
  //Propagate muon and create emited photons
  T.Propagate_Muon();

  //Check if the muon crossed both scintillators (If it did propagate photons)
  if(T.GetDoubleCross())
  {
    T.Propagate_Photons(Simulation.GetN_photons());
    Nphotons_total += Simulation.GetN_photons();
    Nphotons_detected += Simulation.GetN_detected();
    Nphotons_absorbed += Simulation.GetN_absorbed();
    Nphotons_lost += Simulation.GetN_lost();
  }
}*/
