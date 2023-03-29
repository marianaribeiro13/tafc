#include "Generators.h"
#include "Geometry.h"
#include "tools.h"
#include "Particle.h"
#include "Tracker.h"
#include "TApplication.h"

using namespace std;

int main(int argc, char* argv[])
{
  Generator* gen = new Generator();
  double radius = 5.0;
  double height = 1.0;
  double distance = 25.;
  double airgap = 1;
  double althickness = 0.0016;
  double step = 0.001;


  Tracker* T = new Tracker(radius,height,distance,airgap,althickness,step,gen);

  T->Propagate_Muon();
  T->Propagate_Photons();


  TApplication app("app", nullptr, nullptr);
  T->Draw();
  if(argc>1) {app.Run();};

  return 0;
}
