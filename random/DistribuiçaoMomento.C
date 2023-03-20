#include <cmath>
#include <iostream>
#include <vector>
#include "Geometry.h"
#include "Muon.C"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "tools.h"

using namespace std;

int main(int argc, char **argv)
{

    auto f = [](double *x,double *par)
    {
        return pow(cos(x[1]),3)*0.00253*pow(x[0]*cos(x[1]),-(0.2455+1.288*log10(x[0]*cos(x[1]))-0.255*pow(log10(x[0]*cos(x[1])),2)+0.0209*pow(log10(x[0]*cos(x[1])),3) ));
    };
    TF1 *F = new TF1("f",f);
    vector<double> x,y;
    for(int i=0;i<1000;i++)
    {
        vector<double> aux = tools::Random_Distribution_2D(F,0,1000,0,M_PI/2.,F->GetMaximum());
        x.push_back(aux[0]);
        y.push_back(aux[1]);
    }

    auto G = new TGraph(x.size(),x.data(),y.data());
    auto c = new TCanvas("c","c",1000,1000);
    G->Draw("PA");
    c->SaveAs("a.pdf");

    return 0;
}
