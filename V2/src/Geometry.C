#include "Geometry.h"
using namespace std;
///////////////////////////////// Constructor ///////////////////////////////

Geometry::Geometry(double radius, double height, double distance, double airgap, double althickness,int n_SIPMS,double s_size){

    geom = new TGeoManager("telescope", "Telescope geometry");

    Radius = radius;
    Height = height;
    Distance = distance;

    innerradius = radius + airgap;
    outerradius = radius +airgap+althickness;
    Airgap = airgap;
    Thickness = althickness;
    MaxHeight = (distance/2+height+airgap+althickness)*1.2;
    MaxRadius = 2*outerradius;

    n_SIPM = n_SIPMS;
    SIPM_size = s_size;
    SIPM_angle = 2*M_PI / n_SIPM;
    SIPM_alpha = 0.5*SIPM_size/(SIPM_angle*Radius);

    //--- define some materials
    TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0, 0, 0);
    TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98, 13, 2.7);
    TGeoMaterial *matPlastic = new TGeoMaterial("Plastic", 0, 0, 1.023);

    //--- define some media
    TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
    TGeoMedium *Al = new TGeoMedium("Root Material",2, matAl);
    TGeoMedium *Plastic = new TGeoMedium("Root Material2",3, matPlastic);

    //Create world volume
    TGeoVolume* top = geom->MakeTube("TOP", Vacuum, 0, MaxRadius, MaxHeight); // rmin, rmax, mid height
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

    TGeoVolume *sidetube = geom->MakeTube("lateral foil", Al, radius+airgap, radius+airgap+althickness, height/2. +althickness +airgap);
    TGeoVolume *covertube = geom->MakeTube("cover foil", Al, 0, radius+airgap+althickness, althickness/2.);

    TGeoVolume *foil = new TGeoVolumeAssembly("Aluminium foil");
    foil->SetLineColor(kGray);
    foil->AddNode(sidetube, 1, tr3);
    foil->AddNode(covertube, 1, tr4);
    foil->AddNode(covertube, 2, tr5);

    top->AddNode(foil, 1, tr1);
    top->AddNode(foil, 2, tr2);

    geom->CloseGeometry();

}

Geometry::~Geometry()
{
    delete geom;

}

///////////////////////////// Build Muon Telescope //////////////////////////

void Geometry::Build_MuonTelescope(double radius, double height, double distance, double airgap, double althickness,int n_SIPMS,double s_size)
{
    Radius = radius;
    Height = height;
    Distance = distance;

    innerradius = radius + airgap;
    outerradius = radius +airgap+althickness;
    Airgap = airgap;
    Thickness = althickness;
    MaxHeight = (distance/2+height+airgap+althickness)*1.2;
    MaxRadius = 2*outerradius;

    n_SIPM = n_SIPMS;
    SIPM_size = s_size;
    SIPM_angle = 2*M_PI / n_SIPM;
    SIPM_alpha = 0.5*SIPM_size/(SIPM_angle*Radius);

    //--- define some materials
    TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0, 0, 0);
    TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98, 13, 2.7);
    TGeoMaterial *matPlastic = new TGeoMaterial("Plastic", 0, 0, 1.023);

    //--- define some media
    TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
    TGeoMedium *Al = new TGeoMedium("Root Material",2, matAl);
    TGeoMedium *Plastic = new TGeoMedium("Root Material2",3, matPlastic);

    //Create world volume
    TGeoVolume* top = geom->MakeTube("TOP", Vacuum, 0, MaxRadius, MaxHeight); // rmin, rmax, mid height
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

    TGeoVolume *sidetube = geom->MakeTube("lateral foil", Al, radius+airgap, radius+airgap+althickness, height/2. +althickness +airgap);
    TGeoVolume *covertube = geom->MakeTube("cover foil", Al, 0, radius+airgap+althickness, althickness/2.);

    TGeoVolume *foil = new TGeoVolumeAssembly("Aluminium foil");
    foil->SetLineColor(kGray);
    foil->AddNode(sidetube, 1, tr3);
    foil->AddNode(covertube, 1, tr4);
    foil->AddNode(covertube, 2, tr5);

    top->AddNode(foil, 1, tr1);
    top->AddNode(foil, 2, tr2);

    geom->CloseGeometry();

}



bool Geometry::Check_Symmetric_Detector(const double* cpoint)
{
    double h = abs(cpoint[2]);

    if(h>(0.5*(Height+Distance+SIPM_size)) || h<(0.5*(Height+Distance-SIPM_size))){return false;};

    double r = sqrt(cpoint[0]*cpoint[0]+cpoint[1]*cpoint[1]);
    if(abs(r-Radius)>1e-6){return false;};

    double theta = atan(cpoint[1]/cpoint[0])/SIPM_angle;
    double delta = theta - round(theta);
    if(abs(delta) > SIPM_alpha){return false;};

    return true;
}

vector<double> Geometry::GetNormal(double r, double h,vector<double> d,const double* cpoint)
{
  vector<double> aux(3);
  bool horizontal_reflection =false,vertical_reflection=false;


  if(abs(r-Radius) <1e-6 || abs(r-innerradius) <1e-6 || abs(r-outerradius) <1e-6)
  {

    vertical_reflection=true;
    aux[0] = cpoint[0]/r;
    aux[1] = cpoint[1]/r;
    aux[2] = 0;
  }
  if((abs(h-(0.5*Distance+Height))<1e-6) || (abs(h-0.5*Distance)<1e-6) || (abs(h-(Airgap+(0.5*Distance)+Height))<1e-6)  || (abs(h-(-Airgap+(0.5*Distance)))<1e-6) || (abs(h-(Thickness+Airgap+(0.5*Distance)+Height))<1e-6) || (abs(h-(-Thickness-Airgap+(0.5*Distance)))<1e-6)  )
  {
    horizontal_reflection =true;
    aux[0] = 0;
    aux[1] = 0;
    aux[2] = 1;
  }
  if(vertical_reflection && horizontal_reflection)// Corner reflection, photon goes back
  {
    return d;
  }

  if(tools::Angle_Between_Vectors(aux,d)>M_PI/2)
  {
    for(int i=0;i<3;i++){aux[i] = -aux[i];};
  }

  return aux;
}

bool Geometry::VacuumToPlastic(double r,double h)
{
  //If the particle is in vacuum (in the air gap) and near the lateral scintillator boundary
  //or near one of the flat scintillator boundaries
  if(((abs(r-Radius) < 1e-6) || abs(h-(Height/2))<1e-6)){return true;};

  return false;
}

bool Geometry::VacuumToAluminium(double r,double h)
{

  if(((abs(r-innerradius) < 1e-6) || (abs(r-outerradius) <1e-6)
  ||abs(h-(Height/2+Airgap))<1e-6 )){return true;};

  return false;
}
