#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include <iostream>
#include "TCanvas.h"
#include "TGraph2D.h"
#include <cmath>
#include <vector>
#include "TMath.h"

using namespace std;



void Telescope(bool vis = true){


    auto c = new TCanvas("c0","c0",1000,1000);


}


vector<double> generate_vector()
{
    vector<double> v(3);
    v[0] = rand();
    v[1] = rand();
    v[2] = rand();

    double norm = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    v[0] /= norm;
    v[1] /= norm;
    v[2] /= norm;
    return v;
}

int main()
{
    srand(time(0));
    vector<double> x(1000);
    vector<double> y(1000);
    vector<double> z(1000);

    for(int j=0;j<1000;j++)
    {


        cout<< (double)j/100 * 2 * TMath::Pi()<<endl;
        x[j] = 50*cos((double)j/100 * 2 * TMath::Pi());
        y[j] = 25*sin((double)j/100 * 2 * TMath::Pi());
        z[j] = (double) 500-j;
        cout<<x[j]<<" "<<y[j]<<" "<<z[j]<<endl<<endl;

    }

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

    TGeoTranslation *tr1 = new TGeoTranslation(0., 0., 50.5);
    TGeoTranslation *tr2 = new TGeoTranslation(0., 0., -50.5);

    TGeoVolume *scintillator1 = geom->MakeTube("scintillator1", Al, 0, 25., 0.5); // rmin, rmax, mid height
    scintillator1->SetLineColor(kBlue);
    top->AddNode(scintillator1, 1, tr1);
    TGeoVolume *scintillator2 = geom->MakeTube("scintillator2", Al, 0, 25., 0.5); // rmin, rmax, mid height
    scintillator2->SetLineColor(kBlue);
    top->AddNode(scintillator2, 1, tr2);

    geom->CloseGeometry();

    auto G = new TGraph2D(1000,x.data(),y.data(),z.data());
    auto c = new TCanvas("c","c",1600,900);
    G->SetMarkerSize(8);

    top->Draw();
    G->Draw("PSAME");
    c->SaveAs("shi.pdf");


}
