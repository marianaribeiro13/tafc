#include "Modes.h"

void Heatmap_Mode(Tracker* T, DataManager* Data)
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

void EMap_Mode(Tracker* T, DataManager* Data)
{
  vector<double> d = {0,0,-1};
  vector<double> x = {0,0,0};
  x[2] = T->GetHeight()+T->GetDistance()/2 - 1e-12;

  Particle* M0 = new Particle(13,1000,d,x);
  T->Insert_New_Muon(M0);
  T->Propagate_Muon();
  T->Propagate_Photons();
  Data->Fill_Efficiency_Map(T);

  for(int i=0;i<10;i++)
  {
    cout<<i<<endl;
    for(int j=0;j<10;j++)
    {
      x[0] = (double) (i/10) * T->GetRadius() * cos((double) (j/100) * 2 *M_PI);
      x[1] = (double) (i/10) * T->GetRadius() * sin((double) (j/100) * 2 *M_PI);

      Particle* M = new Particle(13,1000,d,x);
      T->Insert_New_Muon(M);
      T->Propagate_Muon();
      T->Propagate_Photons();
      Data->Fill_Efficiency_Map(T);

    }
  }

  Data->Draw_Efficiency_Map("Efficiency_Map.pdf");
  return;
}

void Draw_Mode(Tracker* T, DataManager* Data,int n)
{


  T->Propagate_Muon();
  TApplication app("app", nullptr, nullptr);
  T->Draw(n);
  app.Run();

  return;
}

void Simulation_Mode(Tracker* T,DataManager* Data,char* filename)
{
  for(int i=0;i<1000;i++)
  {
    T->Propagate_Muon();
    if(T->GetDoubleCross())
    {
      T->Propagate_Photons();
    }

    T->Generate_New_Muon();
  }
}
