#include "Geometry.h"


///////////////////////////////// Constructor ///////////////////////////////

Geometry::Geometry(){

    geom = new TGeoManager("telescope", "Telescope geometry");

}


///////////////////////////// Build Muon Telescope //////////////////////////

void Geometry::Build_MuonTelescope(double radius, double height, double distance, double airgap, double althickness){

    //--- define some materials
    TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0, 0, 0);
    TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98, 13, 2.7);
    TGeoMaterial *matPlastic = new TGeoMaterial("Plastic", 0, 0, 1.023);

    //--- define some media
    TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
    TGeoMedium *Al = new TGeoMedium("Root Material",2, matAl);
    TGeoMedium *Plastic = new TGeoMedium("Root Material2",3, matPlastic);
 
    //Create world volume
    TGeoVolume* top = geom->MakeTube("TOP", Vacuum, 0, 200., 100.); // rmin, rmax, mid height
    geom->SetTopVolume(top);

    //Define two transformations to position the scintillators
    TGeoTranslation *tr1 = new TGeoTranslation(0., 0., (double)(distance + height)/2.);
    TGeoTranslation *tr2 = new TGeoTranslation(0., 0., - (double)(distance + height)/2.);

    TGeoVolume *scintillator = geom->MakeTube("scintillator", Plastic, 0, radius, height/2.); // rmin, rmax, mid height
    scintillator->SetLineColor(kOrange+10);
    top->AddNode(scintillator, 1, tr1);
    top->AddNode(scintillator, 2, tr2);

    /////////Create aluminium foil around each scintillator

    //Define transformations to position aluminium foil parts
    TGeoTranslation *tr3 = new TGeoTranslation(0., 0., 0.);
    TGeoTranslation *tr4 = new TGeoTranslation(0., 0., height/2. + airgap + althickness/2.);
    TGeoTranslation *tr5 = new TGeoTranslation(0., 0., -(height/2. + airgap + althickness/2.));

    TGeoVolume *sidetube = geom->MakeTube("lateral foil", Al, radius+airgap, radius+airgap+althickness, height/2. + 2*airgap);
    TGeoVolume *covertube = geom->MakeTube("cover foil", Al, 0, radius+airgap+althickness, althickness/2.);

    TGeoVolume *foil = new TGeoVolumeAssembly("Aluminium foil");
    foil->SetLineColor(kBlack);
    foil->AddNode(sidetube, 1, tr3);
    foil->AddNode(covertube, 1, tr4);
    foil->AddNode(covertube, 2, tr5);

    top->AddNode(foil, 1, tr1);
    top->AddNode(foil, 2, tr2);

    geom->CloseGeometry();

}
