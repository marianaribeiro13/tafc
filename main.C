#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"


int main(){
  TGeoManager *geom = new TGeoManager("telescope", "Telescope geometry");

  //--- define some materials
    TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
    TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98,13,2.7);

    //   //--- define some media
    TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
    TGeoMedium *Al = new TGeoMedium("Root Material",2, matAl);
 
    //--- define the transformations
    TGeoVolume *top = geom->MakeTube("TOP", Vacuum, 0, 50., 100.); // rmin, rmax, mid height
    geom->SetTopVolume(top);

  return 0;
}