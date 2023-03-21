#include "Geometry.h"
#include "Tracking.h"

#include "TApplication.h"

int main(){

  // random generator
  srand(time(0));
  
  // make detector geometry 
  Tracking Track;
  
  /*    
  // generate cosmic muon
  
  Track.AddParticle(...);
  
  // propagate particles
    Track.Propagate();
  */
  
  // display
  TApplication A(nullptr, nullptr);
  Track.Draw();
  A.Run();
  
  return 0;
}
