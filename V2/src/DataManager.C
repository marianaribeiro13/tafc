#include "DataManager.h"

DataManager::DataManager(int n)
{
  x.reserve(n);
  y.reserve(n);
  ZenithAngle.reserve(n);
  N_photons.reserve(n);
  N_absorbed.reserve(n);
  N_lost.reserve(n);
  N_detected.reserve(n);
  DoubleCross.reserve(n);
  Efficiency.reserve(n);
  size = n;
}

DataManager::~DataManager(){}

void DataManager::Extract_Data(Tracker* T,int i)
{
  x[i] = T->GetMuon()->GetStartingPosition()[0];
  y[i] = T->GetMuon()->GetStartingPosition()[1];
  ZenithAngle[i]=tools::Angle_Between_Vectors(T->GetMuon()->GetDirection(),{0,0,-1});
  N_photons[i]=T->GetN_photons();
  N_absorbed[i]=T->GetN_absorbed();
  N_lost[i]=T->GetN_lost();
  N_detected[i]= (double) T->GetN_detected();
  DoubleCross[i] = T->GetDoubleCross();
  Efficiency[i] = (double) N_detected[i]/N_photons[i];
  return;
}

void DataManager::Print_Data(int i)
{
  cout<<"Data from Muon number: "<<i+1<<endl;
  cout<<"Starting Position: "<<x[i]<<" "<<y[i]<<endl;
  cout<<"Zenith Angle: "<<ZenithAngle[i]<<endl;
  cout<<"Total Photons Generated: "<<N_photons[i]<<endl;
  cout<<"Photons Absorbed: "<<N_absorbed[i]<<endl;
  cout<<"Photons Detected: "<<N_detected[i]<<endl;
  cout<<"Photons Lost: "<<N_lost[i]<<endl<<endl;
  return;
}

void DataManager::Draw_Efficiency_Graph()
{
  cout<<x.size()<<endl;
  auto g = new TGraph2D(size,x.data(),y.data(),Efficiency.data());
  auto c = new TCanvas("c","c",1000,1000);
  g->Draw("surf1");
  c->SaveAs("A.pdf");
}
