#include "Modes.h"
#include <thread>
#include <mutex>

using namespace std;


int main(int argc, char* argv[])
{
  auto start = chrono::high_resolution_clock::now();
  double radius = 5.0;
  double height = 1.0;
  double distance = 25.;
  double airgap = .1;
  double althickness = 0.0016;
  double step = 0.001;
  int n_SIPMS = 4;
  double SIPM_size = .6;
  bool flag = false;

  Generator* gen = new Generator();
  DataManager *Data = new DataManager();
  Geometry* Geo = new Geometry(radius,height,distance,airgap,althickness,n_SIPMS,SIPM_size);
  Geo->GetGeoManager()->SetMaxThreads(8);
  // if(argc == 5 && !strcmp(argv[1],"-sim"))
  // {
  //   if(sscanf(argv[2],"%lf",&distance) && sscanf(argv[3],"%d",&n_SIPMS))
  //   {
  //     Tracker* T = new Tracker(radius,height,distance,airgap,althickness,step,gen,n_SIPMS,SIPM_size);
  //     Simulation_Mode(T,Data,argv[4]);
  //     delete T;
  //     goto close;
  //   }
  // }

  Tracker* T = new Tracker(step,gen,Geo);


  if(argc == 2 && !strcmp(argv[1],"-heatmap"))
  {
    Heatmap_Mode(T,Data);
    goto close;
  }
  if(argc==2 && !strcmp(argv[1],"-gEff"))
  {
    GeomEfficiency_Mode(T,Data);
    goto close;
  }
  if(argc == 3 && !strcmp(argv[1],"-heatmap") && !strcmp(argv[2],"1"))
  {
    HeatmapSingle_Mode(T,Data);
    goto close;
  }
  if(argc == 3 && !strcmp(argv[1],"-draw"))
  {
    int n=0;
    if(sscanf(argv[2],"%d",&n))
    {
      Draw_Mode(T,Data,n);
      goto close;
    }
  }

  delete Data;
  delete gen;
  delete T;
  cout<<"Invalid Input"<<endl;
  return -1;

  close:
  delete Data;
  delete T;
  delete gen;

  auto end = chrono::high_resolution_clock::now();
  chrono::duration<float> duration = end-start;
  cout<<"Time Elapsed: "<<duration.count()<<" seconds"<<endl;
  return 0;
}
