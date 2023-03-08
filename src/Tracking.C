#include "Tracking.h"

Tracking::Tracking(Double_t* point,Double_t* dir){

  // Setting both initial point and direction and finding the state:
  TGeoNode* init_node = geom->InitTrack(point, dir);
}

Tracking::~Tracking(){
}

const Double_t* Tracking::GetCurPoint(){
  //Get current point
  const Double_t *cpoint = geom->GetCurrentPoint();
  return cpoint;
}

const Double_t* Tracking::GetCurDir(){
  //Get current direction
  const Double_t *cdir = geom->GetCurrentDirection();
  return cdir;
}

TGeoNode* Tracking::GetCurState(){
  //Get current state (node)
  TGeoNode *current = geom->GetCurrentNode();
  return current;
}

void TVirtualGeoTrack* Tracking::AddTrack(Int_t id, Int_t pdg, TVirtualGeoTrack *parent, TObject* ptrParticle, Double_t xi, Double_t yi, Double_t zi, Double_t ti){
  TVirtualGeoTrack *track = new TVirtualGeoTrack(id, pdgcode, parent, ptrParticle);
  track.AddPoint(xi, yi, zi, ti);
  geom->AddTrack(track)
}

/*
  //Assign current
  track->AddPoint(x,y,z,t);

  //Set arbitrary step
  gGeoManager->SetStep(stepvalue);
  //Execute step and get node(geometrical position) of new position
  TGeoNode *newNode = gGeoManager->Step(kFALSE);
    
  //We probably need to assign the new position to the track
  track->AddPoint(x,y,z,t);
}*/
