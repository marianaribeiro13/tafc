#include "Geometry.h"

Geometry::Geometry(Float_t h, Float_t r, Float_t d) : 
height(h), 
radius(r), 
distance(d) {

    geom = new TGeoManager("telescope", "Telescope geometry");

    //--- define some materials
    TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
    TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98,13,2.7);

    //   //--- define some media
    TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
    TGeoMedium *Al = new TGeoMedium("Root Material",2, matAl);
 
    TGeoVolume *top = geom->MakeTube("TOP", Vacuum, 0, 50., 100.); // rmin, rmax, mid height
    geom->SetTopVolume(top);

    TGeoTranslation *tr1 = new TGeoTranslation(0., 0., float(distance + height)/2.);
    TGeoTranslation *tr2 = new TGeoTranslation(0., 0., float(distance + height)/2.);

    TGeoVolume *scintillator = geom->MakeTube("scintillator", Al, 0, radius, float(height)/2.); // rmin, rmax, mid height
    scintillator->SetLineColor(kBlue);
    top->AddNode(scintillator, 1, tr1);
    top->AddNode(scintillator, 2, tr2);
}