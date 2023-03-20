#include "Tracking.h"
#include "TCanvas.h"
#include "TMath.h"

Tracking::Tracking() : muon_step(0.001), Geometry(25.){

    // Setting both initial point and direction and finding the state:
    // gGeoManager->SetCurrentPoint(x,y,z) + gGeoManager->SetCurrentDirection(nx,ny,nz) + gGeoManager->FindNode()

    //TObject* Muon = new TParticle(pdg,0,0,0,0,0,px,py,pz,E,vx,vy,vz,t);

    //double initial = double(distance + 1.0)/2.;
    std::vector<double> x {0., 0., 0., double(distance + 1.0)/2.};

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

    double Energy = 1.;

    double theta = 2.9 * TMath::Pi()/3.;

    muon = new Muon(x, theta, Energy);

    std::vector<double> d = muon->GetDirection();

    geom->InitTrack(x[1], x[2], x[3], d[0], d[1], d[2]);

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
    track->AddPoint(x[1], x[2], x[3], x[0]);

}

Tracking::~Tracking(){}

void Tracking::Propagation(){

    double alpha = 1.;
    double beta = 1.;

    //double v = 1.;
    double t = 0.;

    while(!(geom->IsOutside())){

        //CrossNextBoundary(muon_step);

        while(geom->IsSameLocation()){
            
            DefinedStep(muon_step);

            double v = muon->GetVelocity();

            t += double(muon_step)/v;

            double E = muon->GetEnergy();

            double dE = (-1)*alpha + (-1)*beta * E;

            muon->ChangeEnergy(E-dE);
            muon->ChangeMomentum();
            muon->ChangeVelocity();

            const double *cpoint = geom->GetCurrentPoint();
    
            //Assign the new position to the track
            track->AddPoint(*cpoint, *(cpoint+1), *(cpoint+2),t);
        }

        double vacuum_step = CrossNextBoundary();

        t += double(vacuum_step)/(muon->GetVelocity());

        const double *cpoint = geom->GetCurrentPoint();

        track->AddPoint(*cpoint, *(cpoint+1), *(cpoint+2),t);

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

void Tracking::Drawing(){

    TCanvas *c1 = new TCanvas("c1","c1");
    c1->cd();
    //top->Draw("ogle");
    top->Draw();
    geom->DrawTracks();
    //geom->AnimateTracks(0, 1000, 200, "/* /G /S");
    c1->SaveAs("Drawing.pdf");
}
