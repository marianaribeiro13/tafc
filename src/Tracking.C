#include "Tracking.h"
#include "TCanvas.h"

Tracking::Tracking() : muon_step(0.02), Geometry(25.){

    // Setting both initial point and direction and finding the state:
    // gGeoManager->SetCurrentPoint(x,y,z) + gGeoManager->SetCurrentDirection(nx,ny,nz) + gGeoManager->FindNode()

    Int_t pdg = 13;
    Double_t c = 3.0 * pow(10,8);
    Double_t x = 0;
    Double_t y = 0;
    Double_t z = double(distance + 1.0)/2.;
    Double_t px = 1.;
    Double_t py = 1.;
    Double_t pz = 1.;
    Double_t E = 1.;
    Double_t vx = 1.;
    Double_t vy = 1.;
    Double_t vz = -1.;
    Double_t t = 0.;
    Double_t v = sqrt(vx*vx + vy*vy + vz*vz);
    Double_t nx = double(vx)/v;
    Double_t ny = double(vy)/v;
    Double_t nz = double(vz)/v;
    Double_t gamma = 1./sqrt(1-v*v/c*c);

    //TObject* Muon = new TParticle(pdg,0,0,0,0,0,px,py,pz,E,vx,vy,vz,t);

    //TParticle* Muon = new TParticle();

    TObject* Muon = new TObject();

    geom->InitTrack(x, y, z, nx, ny, nz);

    //Const Double_t *firstpoint = geom->GetCurrentPoint();

    /*//Get current point
    Const Double_t *cpoint = geom->GetCurrentPoint();
    //Get current direction
    Const Double_t *cdir = geom->GetCurrentDirection();
    //Get current state (node)
    TGeoNode *current = geom->GetCurrentNode();*/

    //Create Track associated to particle
    Int_t track_index = geom->AddTrack(0,pdg,Muon); //track id, particle pdg, pointer to particle

    //Set the created track as the current track
    geom->SetCurrentTrack(track_index);
    //Get pointer to current track
    track = geom->GetCurrentTrack();
    //Assign current 
    track->AddPoint(x, y, z, t);

}

Tracking::~Tracking(){}

void Tracking::Propagation(){

    //Double_t v = sqrt(Muon->Vx()*Muon->Vx()+Muon->Vy()*Muon->Vy()+Muon->Vz()*Muon->Vz());

    Double_t v = 1;

    Double_t t = 0;

    while(!(geom->IsOutside())){

        CrossNextBoundary(muon_step);

        const Double_t *cpoint = geom->GetCurrentPoint();

        t += double(muon_step)/v;
    
        //Assign the new position to the track
        track->AddPoint(*cpoint, *(cpoint+1), *(cpoint+2),t);
    }
}

void Tracking::CrossNextBoundary(Double_t pstep){

    //Sets the possible step size
    geom->FindNextBoundary(pstep);

    //Gets the biggest step possible in any direction that assures that no boundary will be crossed
    //Double_t safety = geom->GetSafeDistance();

    //Get step value
    Double_t snext = geom->GetStep();
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

    //return newNode;
}

TGeoNode* Tracking::DefinedStep(Double_t stepvalue){

    //Set arbitrary step
    geom->SetStep(stepvalue);
    //Execute step and get node(geometrical position) of new position
    TGeoNode *newNode = geom->Step(kFALSE);

    return newNode;
}

void Tracking::Drawing(){

    TCanvas *c1 = new TCanvas("c1","c1");
    c1->cd();
    //top->Draw("ogle");
    top->Draw();
    geom->DrawTracks();
    //geom->AnimateTracks(0, 1000, 200, "/G /S");
    c1->SaveAs("please.pdf");
}
