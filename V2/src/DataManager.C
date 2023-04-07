#include "DataManager.h"

DataManager::DataManager()
{
  Heatmap = new TH2D("h","",20,-M_PI,M_PI,20,0,1);
  Efficiency_Map = new TH2D("h","",20,-5,5,20,-5,5);

}

DataManager::~DataManager()
{
  delete Heatmap;
  delete Efficiency_Map;
}

void DataManager::Fill_Efficiency_Map(Tracker* T)
{
  double efficiency = (double)T->GetN_detected()/T->GetN_photons();
  Efficiency_Map->Fill(T->GetMuon()->GetStartingPosition()[0],T->GetMuon()->GetStartingPosition()[1],efficiency);
}


void DataManager::Draw_Efficiency_Map(string name)
{
  auto c = new TCanvas("c","c",1600,1000);
  Efficiency_Map->Draw("COLZPOL");
  c->SaveAs(name.c_str());
  delete c;
}

void DataManager::Fill_Heatmap(Tracker* T)
{
  if(T->GetPhotonFlag()){return;};
  T->Fill_Heatmap(Heatmap);

  return;
}

void DataManager::Draw_Heatmap(string name)
{
  auto c = new TCanvas("c","c",1600,1000);
  Heatmap->Draw("colz");
  c->SaveAs(name.c_str());
  delete c;

  return;
}
