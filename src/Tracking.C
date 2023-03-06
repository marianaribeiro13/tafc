#include "Tracking.h"

Tracking::Tracking(){

    // Setting both initial point and direction and finding the state:
    geom->InitTrack(Double_t *point[3],Double_t *dir[3]);

    //Get current point
    Const Double_t *cpoint = gGeoManager->GetCurrentPoint();
    //Get current direction
    Const Double_t *cdir = gGeoManager->GetCurrentDirection();
    //Get current state (node)
    TGeoNode *current = gGeoManager->GetCurrentNode();


    //Create Track associated to particle
    Int_t track_index = gGeoManager->AddTrack(id,pdg,ptrParticle);

    //Set the created track as the current track
    gGeoManager->SetCurrentTrack(track_index);
    //Get pointer to current track
    TVirtualGeoTrack *track = gGeoManager->GetCurrentTrack();
    //Assign current 
    track->AddPoint(x,y,z,t);

    //Set arbitrary step
    gGeoManager->SetStep(stepvalue);
    //Execute step and get node(geometrical position) of new position
    TGeoNode *newNode = gGeoManager->Step(kFALSE);
    
    //We probably need to assign the new position to the track
    track->AddPoint(x,y,z,t);
}