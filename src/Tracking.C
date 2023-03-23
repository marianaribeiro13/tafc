#include "Tracking.h"
#include "TCanvas.h"
#include "TMath.h"
#include <cmath>

using namespace std;
//NATURAL
Double_t alpha=0.0072974;
Double_t I1=64.70713358e-6; //MeV
Double_t n=2.5587248e18; //MeV
Int_t Z=-1;
Double_t mass_muon=105.6583755; //MeV
Double_t mass_electron=0.51099895; //MeV
Double_t dE=0;
Double_t dx=1.e-4;

//S.I.
Double_t c = 299792458;
Double_t mp = 1.672621637e-27;
Double_t me = 9.1093821499999992e-31;
Double_t qe = 1.602176487e-19;
Double_t na = 6.02214179e23;
Double_t eps0 = 8.854187817e-12;
Double_t n_density=3.33e29; // per m^3
Double_t I = 1.03660828e-17;

Double_t dEdx_nat(Double_t v) //Mev/cm
{
  return (4*M_PI*n*Z*Z*alpha*alpha*(log((2*mass_electron*v*v)/(I1*(1-v*v)))-v*v))/(mass_electron*v*v);
}

Double_t dEdx_func(Double_t beta) //S.I.
{
  Double dE_SI;
  dE_SI = (qe*qe*qe*qe*n_density*Z*Z*(log((2*me*c*c*beta*beta)/(I*(1-beta*beta)))-beta*beta))/(4*M_PI*eps0*eps0*me*c*c*beta*beta);
  dE_SI=dE_SI/(1.602e-13);
  return dE_SI;
}  

Tracking::Tracking() : muon_step(0.001), Geometry(25.){

  Generator* G= Generator();

    // Setting both initial point and direction and finding the state:
    // gGeoManager->SetCurrentPoint(x,y,z) + gGeoManager->SetCurrentDirection(nx,ny,nz) + gGeoManager->FindNode()

    //TObject* Muon = new TParticle(pdg,0,0,0,0,0,px,py,pz,E,vx,vy,vz,t);

    //double initial = double(distance + 1.0)/2.;
    std::vector<double> x = G->Generate_Position(5.0,distance,1.0);

    /*auto spec = [](double* p,double* par) {

	    return p[0]*par[0]*pow(p[0],(-1)*(par[1] + par[2]*log10(p[0]) + par[3]*log10(p[0])*log10(p[0]) + par[4]*log10(p[0])*log10(p[0])*log10(p[0])));
    };

    TF1 MuonSpectrum("MuonSpectrum", spec, 1, 20, 0); //npar=0,ndim=1(default)

    MuonSpectrum.SetParameter(0.00253,0.2455,1.288,-0.2555,0.0209)

    auto gaussian = [](double* p,double* par) {

	    return (1./sqrt(2*TMath::Pi()*par[1]*par[1])) * exp((-1)*((p-par[0])*(p-par[0]))/(2*par[1]*par[1]));
    };

    TF1 Gaussian("Gaussian", gaussian, 1, 20, 0); //npar=0,ndim=1(default)

    Gaussian.SetParameter(0,1);

    auto invgaussian = [](double* p,double* par) {

	    return (1./sqrt(2*TMath::Pi()*par[1]*par[1])) * exp((-1)*((p-par[0])*(p-par[0]))/(2*par[1]*par[1]));
    };

    TF1 InvGaussian("InvGaussian", invgaussian, 1, 20, 0); //npar=0,ndim=1(default)

    Gaussian.SetParameter(0,1);


    IntegratorMC I(MuonSpectrum);*/
    
    //double momentum = ;
    auto f = [](double *x,double *par)
    {
        return pow(cos(x[1]),3)*0.00253*pow(x[0]*cos(x[1]),-(0.2455+1.288*log10(x[0]*cos(x[1]))-0.255*pow(log10(x[0]*cos(x[1])),2)+0.0209*pow(log10(x[0]*cos(x[1])),3) ));
    };

    TF1 *F = new TF1("f",f);

    std::vector<double> aux = G->Random_Distribution_2D(F,1,2000,0,M_PI/2,F->GetMaximum());

    double Momentum = aux[0];

    double theta = aux[1];

    muon = new Muon(x, theta, Momentum);

    std::vector<double> d = muon->GetDirection();

    geom->InitTrack(x[0], x[1], x[2], d[0], d[1], d[2]);

    //Const Double_t *firstpoint = geom->GetCurrentPoint();*/

    /*//Get current point
    Const Double_t *cpoint = geom->GetCurrentPoint();
    //Get current direction
    Const Double_t *cdir = geom->GetCurrentDirection();
    //Get current state (node)
    TGeoNode *current = geom->GetCurrentNode();*/

    //Create Track associated to particle
    int track_index = geom->AddTrack(0, muon->GetPDG(), muon); //track id, particle pdg, pointer to particle

    //Set the created track as the current track
    geom->SetCurrentTrack(track_index);
    //Get pointer to current track
    track = geom->GetCurrentTrack();
    //Assign current 
    track->AddPoint(x[0], x[1], x[2], 0);

}

Tracking::~Tracking(){}

void Tracking::Propagate(){

    //double alpha = 1.;
    //double beta = 1.;

    //double v = 1.;
    double t = 0.;

    while(!(geom->IsOutside())){

        //CrossNextBoundary(muon_step);

        double air_step = CrossNextBoundary();

        t += double(air_step)/(muon->GetVelocity());

        const double *cpoint = geom->GetCurrentPoint();

        track->AddPoint(*cpoint, *(cpoint+1), *(cpoint+2),t);

        //std::cout << "Enter Scintillator" << std::endl;

        DefinedStep(0.000001);

        Double_t x=0;
        Double_t dE=0;
        //Double_t dx=1.e-4;

        //std::cout << geom->IsSameLocation() << std::endl;

        while(geom->IsSameLocation()){
            
            DefinedStep(muon_step);

            double v = muon->GetVelocity();

            t += double(muon_step)/v;

            double E = muon->GetEnergy();

            cout << E << endl;

            dE=dEdx_func(v)*muon_step;

            //double dE = ((-1)*alpha + (-1)*beta * E)*muon_step;

            muon->ChangeEnergy(E-dE);
            muon->ChangeMomentum();
            muon->ChangeVelocity();

            const double *cpoint = geom->GetCurrentPoint();

            //Assign the new position to the track
            track->AddPoint(*cpoint, *(cpoint+1), *(cpoint+2),t);
            //}

        }

        //std::cout << "Air" << std::endl;

    }
}

double Tracking::CrossNextBoundary(){

    //Sets the possible step size
    geom->FindNextBoundary();

    //Gets the biggest step possible in any direction that assures that no boundary will be crossed
    //Double_t safety = geom->GetSafeDistance();

    //Get step value
    double snext = geom->GetStep();
    //The geometrical step is taken
    TGeoNode *newNode = geom->Step();
    // The step=snext+epsil is made

    /*Bool_t hasCrossed = gGeoManager->IsEntering();
    // Is the boundary crossed or not?
    Bool_t isOnBoundary = gGeoManager->IsOnBoundary(); // The proposed
    // geometrically limited step to be made was smaller
    // than epsil value.
    Bool_t isOutside = gGeoManager->IsOutside();
    //Did we exit geometry ?*/

    return snext;
}

void Tracking::DefinedStep(double stepvalue){

    //Set arbitrary step
    geom->SetStep(stepvalue);
    //Execute step and get node(geometrical position) of new position
    TGeoNode *newNode = geom->Step(kFALSE);

}

void Tracking::Draw(){

    TCanvas *c1 = new TCanvas("c1","c1");
    c1->cd();
    //top->Draw("ogle");
    top->Draw();
    geom->DrawTracks();
    //geom->AnimateTracks(0, 1000, 200, "/* /G /S");
    c1->SaveAs("Drawing.pdf");
}
