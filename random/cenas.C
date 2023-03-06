#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include <iostream>
#include "TCanvas.h"
#include <cmath>
#include "TF1.h"
#include "TGraph.h"
#include "TAxis.h"
#include <vector>
#include "TMath.h"
#include "Math/IntegratorMultiDim.h"

using namespace std;



int main()
{
    auto II = [](double *x,double *par)
    {
        return cos(x[1])*cos(x[1])*cos(x[1])*0.00253*pow(x[0],-(0.2455 + 1.288*log10(x[0]*cos(x[1])) -0.2555*log10(x[0]*cos(x[1]))*log10(x[0]*cos(x[1])) + 0.0209*log10(x[0]*cos(x[1]))*log10(x[0]*cos(x[1]))*log10(x[0]*cos(x[1]))));;
    };
    auto I = new TF1("Momentum distribution",II,0,20,0);
    auto F = new ROOT::Math::IntegratorMultiDim();
    auto c = new TCanvas("c","c",1000,1000);
    I->Draw();
    c->SaveAs("shi.pdf");

    cout<<I->Eval(5,0)<<endl;
    return 0;
}
