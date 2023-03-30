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
  double radius = 5.0;
  double height = 1.0;
  double distance = 25.;
  double airgap = 1;
  double althickness = 0.0016;
  double step = 0.001;

  if(argc > 3 || (strncmp(argv[1],"-debug",100) && strncmp(argv[1],"-draw",100)) )
  {
    cout<<"Invalid Input"<<endl;
    return 1;
  }

  if(argc == 3 && !strncmp(argv[1],"-debug",100)) //Debug mode
  {
    int seed=0;
    if(sscanf(argv[2],"%d",&seed))
    {
      Generator* gen = new Generator(seed);

      Tracker* T = new Tracker(radius,height,distance,airgap,althickness,step,gen);

      T->Debug();
      cout<<"Generator Seed: "<<seed<<endl;
    }
  }

  if(argc == 3 && !strncmp(argv[1],"-draw",100)) //Draw mode
  {
    int n=0;
    if(sscanf(argv[2],"%d",&n))
    {
      Generator* gen = new Generator();

      Tracker* T = new Tracker(radius,height,distance,airgap,althickness,step,gen);

      T->Propagate_Muon();
      T->Propagate_Photons_DrawMode(n);

      TApplication app("app", nullptr, nullptr);
      T->Draw();
      app.Run();
    }

  }

  if(argc ==1) //Main program
  {
    Generator* gen = new Generator();

    Tracker* T = new Tracker(radius,height,distance,airgap,althickness,step,gen);

    T->Propagate_Muon();
    T->Propagate_Photons();
  }

  return 0;
}
