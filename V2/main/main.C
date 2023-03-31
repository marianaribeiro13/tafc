#include "Generators.h"
#include "Geometry.h"
#include "tools.h"
#include "Particle.h"
#include "Tracker.h"
#include "TApplication.h"
#include <string>

using namespace std;

int main(int argc, char* argv[])
{
  int t0 = time(0),t1=0;
  double radius = 5.0;
  double height = 1.0;
  double distance = 25.;
  double airgap = .1;
  double althickness = 0.0016;
  double step = 0.001;
  int n_SIPMS = 4;

  if(argc == 1) //Simulates a single muon
  {
    Generator* gen = new Generator();

    Tracker* T = new Tracker(radius,height,distance,airgap,althickness,step,gen,n_SIPMS);
    T->Propagate_Muon();
    T->Propagate_Photons();

    cout<<endl<<"Total Photons Generated: "<<T->GetN_photons()<<endl;
    cout<<"Photons Absorbed: "<<T->GetN_absorbed()<<endl;
    cout<<"Photons Detected: "<<T->GetN_detected()<<endl;
    cout<<"Photons Lost: "<<T->GetN_lost()<<endl;

    t1 = time(0);
    cout<<"Time elapsed: "<<t1-t0<<" seconds"<<endl;
    cout<<"Muons Simulated: "<<1<<endl;
    return 0;
  }

  if(argc == 2) //Simulates n muons
  {
    int n_muons=0;
    if(sscanf(argv[1],"%d",&n_muons))
    {
      int N_photons=0,N_absorbed=0,N_detected=0,N_lost=0;
      Generator* gen = new Generator();
      for(int i=0;i<n_muons;i++)
      {
        cout<<endl<<"Muon: "<<i<<endl<<endl;

        Tracker* T = new Tracker(radius,height,distance,airgap,althickness,step,gen,n_SIPMS);
        T->Propagate_Muon();
        T->Propagate_Photons();
        N_photons+=T->GetN_photons();
        N_absorbed+=T->GetN_absorbed();
        N_detected+=T->GetN_detected();
        N_lost+=T->GetN_lost();
      }

      cout<<endl<<"Total Photons Generated: "<<N_photons<<endl;
      cout<<"Photons Absorbed: "<<N_absorbed<<endl;
      cout<<"Photons Detected: "<<N_detected<<endl;
      cout<<"Photons Lost: "<<N_lost<<endl;

      t1 = time(0);
      cout<<"Time elapsed: "<<t1-t0<<" seconds"<<endl;
      cout<<"Muons Simulated: "<<n_muons<<endl;
      return 0;
    }
  }

  if(argc == 3 && !strncmp(argv[1],"-debug",100)) //Debug mode - Prints information about every photon
  {
    int seed=0;
    if(sscanf(argv[2],"%d",&seed))
    {
      Generator* gen = new Generator(seed);

      Tracker* T = new Tracker(radius,height,distance,airgap,althickness,step,gen,n_SIPMS);

      T->Debug();

      TApplication app("app", nullptr, nullptr);
      T->Draw();
      app.Run();

      cout<<"Generator Seed: "<<seed<<endl;
      return 0;
    }

  }

  if(argc == 3 && !strncmp(argv[1],"-draw",100)) //Draw mode - Single Muon and draw n photons
  {
    int n=0;
    if(sscanf(argv[2],"%d",&n))
    {
      Generator* gen = new Generator();

      Tracker* T = new Tracker(radius,height,distance,airgap,althickness,step,gen,n_SIPMS);

      T->Propagate_Muon();
      T->Propagate_Photons_DrawMode(n);

      TApplication app("app", nullptr, nullptr);
      T->Draw();
      app.Run();

      cout<<endl<<"Total Photons Generated: "<<T->GetN_photons()<<endl;
      cout<<"Photons Absorbed: "<<T->GetN_absorbed()<<endl;
      cout<<"Photons Detected: "<<T->GetN_detected()<<endl;
      cout<<"Photons Lost: "<<T->GetN_lost()<<endl;

      return 0;
    }

  }

  if(argc == 3 && !strncmp(argv[1],"-seed",100)) //Seed mode - Single Muon with a given seed
  {
    int seed=0;
    if(sscanf(argv[2],"%d",&seed))
    {
      Generator* gen = new Generator(seed);

      Tracker* T = new Tracker(radius,height,distance,airgap,althickness,step,gen,n_SIPMS);

      T->Propagate_Muon();
      T->Propagate_Photons();

      cout<<endl<<"Total Photons Generated: "<<T->GetN_photons()<<endl;
      cout<<"Photons Absorbed: "<<T->GetN_absorbed()<<endl;
      cout<<"Photons Detected: "<<T->GetN_detected()<<endl;
      cout<<"Photons Lost: "<<T->GetN_lost()<<endl;

      return 0;
    }

  }

  cout<<"Invalid Input"<<endl;
  return 1;
}
