#include "Geometry.h"

Geometry::Geometry(double d) :  
distance(d) {

    geom = new TGeoManager("telescope", "Telescope geometry");

    //--- define some materials
    TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
    TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98,13,2.7);

    //   //--- define some media
    TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
    TGeoMedium *Al = new TGeoMedium("Root Material",2, matAl);
 
    top = geom->MakeTube("TOP", Vacuum, 0, 200., 100.); // rmin, rmax, mid height
    geom->SetTopVolume(top);

    TGeoTranslation *tr1 = new TGeoTranslation(0., 0., double(distance + 1.0)/2.);
    TGeoTranslation *tr2 = new TGeoTranslation(0., 0., - double(distance + 1.0)/2.);

    TGeoVolume *scintillator = geom->MakeTube("scintillator", Al, 0, 5.0, 0.5); // rmin, rmax, mid height
    scintillator->SetLineColor(kBlue);
    top->AddNode(scintillator, 1, tr1);
    top->AddNode(scintillator, 2, tr2);

    geom->CloseGeometry();
}

Geometry::~Geometry(){}