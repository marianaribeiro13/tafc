#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"


int main(){
  new TGeoManager("world", "A simples geometry.");
  TGeoMaterial *mat = new TGeoMaterial("Vacuum",0,0,0);
  TGeoMedium *med = new TGeoMedium("Vacuum",1,mat);
  TGeoVolume *top=gGeoManager->MakeTube("Top",med,10.,10.,10.);
  gGeoManager->SetTopVolume(top);
  gGeoManager->CloseGeometry();
  top->SetLineColor(kMagenta);
  gGeoManager->SetTopVisible();
  top->Draw();
  return 0;
}
