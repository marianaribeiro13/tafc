#include "Tracker.h"

///////////////////////////// Constructor //////////////////////////////

Tracker::Tracker(double radius, double height, double distance, double airgap, double althickness, double step, Generator* g, int n_SIPMS,double SIPM_size)
: Geometry(), stepvalue(step), N_photons(0),N_absorbed(0),N_detected(0),N_lost(0),DoubleCross(false)
{
  generator = g;

  //Build the Telescope Geometry
  Build_MuonTelescope(radius, height, distance, airgap, althickness,n_SIPMS,SIPM_size);

  //Generate Cosmic Muon at incident scintillator plane
  Muon = generator->Generate_CosmicMuon(generator->Generate_Position(Distance,Height,Radius));

  //BetheBloch lambda function
  auto f = [](double *x,double *par)
  {
    int Z=-1;
    double c = 299792458;
    double mp = 1.672621637e-27;
    double me = 9.1093821499999992e-31;
    double qe = 1.602176487e-19;
    double na = 6.02214179e23;
    double eps0 = 8.854187817e-12;
    double n_density=3.33e29; // per m^3
    double I = 1.03660828e-17;
    return ((qe*qe*qe*qe*n_density*Z*Z*(log((2*me*c*c*x[0]*x[0])/(I*(1-x[0]*x[0])))-x[0]*x[0]))/(4*M_PI*eps0*eps0*me*c*c*x[0]*x[0]))/(1.602e-13);
  };
  BetheBloch = new TF1("f",f); //Create BetheBloch TF1

  //Add muon track to the Navigator and assign it to main_track atribute
  geom->AddTrack(0,Muon->GetPDG(),Muon);
  main_track = geom->GetTrack(0);

  //Add starting point to the track (initial position of the muon) and set main_track as the current track
  main_track->AddPoint(Muon->GetStartingPosition()[0],Muon->GetStartingPosition()[1],Muon->GetStartingPosition()[2],0);
  geom->SetCurrentTrack(0);

  //Setting both initial point and direction and finding the state
  geom->InitTrack(geom->GetCurrentTrack()->GetFirstPoint(), Muon->GetDirection().data());

  //Get pointer to current position (Whenever we update the position cpoint is updated,
  //since it is equal to the geom pointer to the current position)
  cpoint = geom->GetCurrentPoint();
}

///////////////////////////////// Destructor ///////////////////////////////////

Tracker::~Tracker()
{

  delete Muon;
}

//Update energy and momentum of the particle after energy loss by BetheBloch equation
double Tracker::Update_Energy(double step)
{
  double dE = BetheBloch->Eval(Muon->GetVelocity()) * (step/100); //step converted from cm to m
  Muon->ChangeEnergy(Muon->GetEnergy()-dE); //Update energy
  Muon->ChangeMomentum(Muon->CalculateMomentum(Muon->GetEnergy())); //Update momentum
  return dE;
}

//Checks whether the next position (if we make the defined step value in the current particle direction)
//is in the same node/volume of the current position
bool Tracker::CheckSameLocation()
{
  double aux[3];
  for(int i=0;i<3;i++)
  {
    aux[i] = cpoint[i] + stepvalue*Muon->GetDirection()[i];
  }
  return geom->IsSameLocation(aux[0],aux[1],aux[2]);
}

//Calculate probability of reflection at boundary
double Tracker::FresnelLaw(double thetai, double n1, double n2)
{

    // Reflection probability for s-polarized light
    double Rs = abs((n1*cos(thetai)-n2*sqrt(1-(n1*sin(thetai)/n2)*(n1*sin(thetai)/n2)))/
                    (n1*cos(thetai)+n2*sqrt(1-(n1*sin(thetai)/n2)*(n1*sin(thetai)/n2))));

    // Reflection probability for p-polarized light
    double Rp = abs((n1*sqrt(1-(n1*sin(thetai)/n2)*(n1*sin(thetai)/n2))-n2*cos(thetai)))/
                    (n1*sqrt(1-(n1*sin(thetai)/n2)*(n1*sin(thetai)/n2))+n2*cos(thetai));

    return 0.5*(Rs+Rp);
}

//Check if light is reflected or transmitted and get new light direction
bool Tracker::CheckReflection(double thetai,double n1,double n2){

  //Total reflection
  if(n1>n2 && thetai > asin(n2/n1))
  {
    return true;
  }

  double Reff = FresnelLaw(thetai, n1, n2);

  if(generator->Uniform(0,1) < Reff) {
    //Photon reflected
    return true;

  } else {
    //Photon transmited
    return false;
  }
}

//Get normal vector to next crossing surface
vector<double> Tracker::GetNormal()
{
  double r = sqrt(cpoint[0]*cpoint[0]+cpoint[1]*cpoint[1]);
  double h = abs(cpoint[2]);
  vector<double> aux(3);
  bool horizontal_reflection =false,vertical_reflection=false;
  vector<double> d = {geom->GetCurrentDirection()[0],geom->GetCurrentDirection()[1],geom->GetCurrentDirection()[2]};

  if(abs(r-Radius) <1e-6 || abs(r-innerradius) <1e-6 || abs(r-outerradius) <1e-6) //If this condition is true, the photon is being reflected by the lateral disk wall
  {

    vertical_reflection=true;
    aux[0] = cpoint[0]/r;
    aux[1] = cpoint[1]/r;
    aux[2] = 0;
  }
  //If this condition is true, the photon is being reflected off of a horizontal scintillator boundary
  if((abs(h-(0.5*Distance+Height))<1e-6) || (abs(h-0.5*Distance)<1e-6) || (abs(h-(Airgap+(0.5*Distance)+Height))<1e-6)  || (abs(h-(-Airgap+(0.5*Distance)))<1e-6) || (abs(h-(Thickness+Airgap+(0.5*Distance)+Height))<1e-6) || (abs(h-(-Thickness-Airgap+(0.5*Distance)))<1e-6)  )
  {
    horizontal_reflection =true;
    aux[0] = 0;
    aux[1] = 0;
    aux[2] = 1;
  }
  if(vertical_reflection && horizontal_reflection)// If the photon happens to be exactly at the corner of the scintillator it will be reflected both horizontally and vertically
  {
    return d;
  }

  if(tools::Angle_Between_Vectors(aux,d)>M_PI/2) //We define the normal vector to have a positive dot product with the direction vector
  {
    for(int i=0;i<3;i++){aux[i] = -aux[i];};
  }

  return aux;
}

//Check if the particle is in vacuum near the scintillator boundary
bool Tracker::VacuumToPlastic(double r)
{
  //If the particle is in vacuum (in the air gap) and near the lateral scintillator boundary
  //or near one of the flat scintillator boundaries
  if(CheckDensity()==0 && ((abs(r-Radius) < 1e-6) || (abs(abs(cpoint[2])-(0.5*Distance+Height))<1e-6)
  || (abs(abs(cpoint[2])-(0.5*Distance))<1e-6))){return true;};

  return false;
}

//Check if the particle is in vacuum near aluminium foil
bool Tracker::VacuumToAluminium(double r)
{
  //If the particle is in vacuum (either in the air gap or outside the telescope) and near any aluminium foil boundary
  if(CheckDensity()==0 && ((abs(r-innerradius) < 1e-6) || (abs(r-outerradius) <1e-6)
  || (abs(abs(cpoint[2])-(Airgap+0.5*Distance+Height))<1e-6) || (abs(abs(cpoint[2])-(-Airgap+(0.5*Distance)))<1e-6)
  || (abs(abs(cpoint[2])-(Thickness+Airgap+0.5*Distance+Height))<1e-6) || (abs(abs(cpoint[2])-(-Thickness-Airgap+(0.5*Distance)))<1e-6) ) )
  {return true;};

  return false;
}

//Check if photon is in the detector region and if it is detected according to the detector efficiency
bool Tracker::DetectionCheck(const double *cpoint,double e)
{
  if(Check_Symmetric_Detector(cpoint) && generator->Uniform(0,1)<generator->GetDetector_Efficiency()->Eval(e))
  {
    return true;
  }
  return false;
}


//////////Muon Propagators//////////
//////////Muon Propagators//////////
//////////Muon Propagators//////////
//////////Muon Propagators//////////
//////////Muon Propagators//////////

void Tracker::Propagate_Muon()
{
  //Variable that stores the number of times the particle crosses a scintillator
  int scintillator_cross=0;

  while(!geom->IsOutside()) //While the particle is inside the defined top volume
  {
    if(CheckDensity()==0) //Particle is in vacuum
    {
      Muon_Vacuum_Step();
    }
    if(CheckDensity()==1.023) //Particle is in a scintillator
    {
      scintillator_cross++; //Plus one scintillator crossed
      Muon_Scintillator_Step();
    }
    if(CheckDensity()==2.7) //Particle is in the aluminium foil
    {
      Muon_Aluminium_Step();
    }
    if(scintillator_cross == 2){DoubleCross=true;}; //The particle crosses both scintillators

  }
  return;
}

//Make vacuum step
void Tracker::Muon_Vacuum_Step()
{
  geom->FindNextBoundaryAndStep(); //Make step to the next boundary and cross it (updates current position of the navigator)
  Muon->IncreaseTime(geom->GetStep()/(Muon->GetVelocity()*2.998e10)); //Update time
  Muon->ChangePosition(cpoint); //Update muon object position
  geom->GetCurrentTrack()->AddPoint(cpoint[0],cpoint[1],cpoint[2],Muon->GetTime()); //Add new point to the current track
  return;
}

//Make scintillator steps until the particle leaves the scintillator
void Tracker::Muon_Scintillator_Step()
{
  //Set next step to the defined arbitrary step value (defined by the constructer)
  geom->SetStep(stepvalue);
  while(CheckSameLocation())
  {
    //Execute step defined by SetStep (flag = kFALSE means it is an arbitrary step, not limited by geometrical reasons)
    geom->Step(kFALSE); //Updates the current position of the navigator
    Muon->IncreaseTime(stepvalue/(Muon->GetVelocity()*2.998e10)); //Update time
    Muon->ChangePosition(cpoint); //Update muon object position
    geom->GetCurrentTrack()->AddPoint(cpoint[0],cpoint[1],cpoint[2],Muon->GetTime()); //Add new point to the current track

    //Get number of photons from Poisson distribution (usually approximatly 1 photon per 100 eV for common scintillators)
    int n = generator->Generate_Photon_Number(10000*Update_Energy(stepvalue)); //Energy is converted to MeV
    N_photons+=n; //Update total number of emited photons

    for(int i=0;i<n;i++)
    {
      //Add generated photon according to Scintillator spectrum to the vector of photons (to propagate later)
      Photons.push_back(generator->Generate_Photon(Muon->GetPosition()));
      Photons[i]->IncreaseTime(Muon->GetTime()); //Assign current time to the photon
    }
  }

  geom->FindNextBoundaryAndStep(stepvalue); //Make step to the next boundary and cross it (step is less than the defined step)
  Muon->IncreaseTime(geom->GetStep()/(Muon->GetVelocity()*2.998e10)); //Update time
  Muon->ChangePosition(cpoint); //Update muon object position
  geom->GetCurrentTrack()->AddPoint(cpoint[0],cpoint[1],cpoint[2],Muon->GetTime()); //Add new point to the current track

  //Get number of photons from Poisson distribution (usually approximatly 1 photon per 100 eV for common scintillators)
  int n = generator->Generate_Photon_Number(10000*Update_Energy(geom->GetStep()));
  N_photons+=n; //Update total number of emited photons

  for(int i=0;i<n;i++)
  {
    //Add generated photon according to Scintillator spectrum to the vector of photons (to propagate later)
    Photons.push_back(generator->Generate_Photon(Muon->GetPosition()));
    Photons[i]->IncreaseTime(Muon->GetTime()); //Assign current time to the photon
  }

  return;
}

//Make aluminium step
void Tracker::Muon_Aluminium_Step()
{
  geom->FindNextBoundaryAndStep(); //Make step to the next boundary and cross it (updates current position of the navigator)
  Muon->IncreaseTime(geom->GetStep()/(Muon->GetVelocity()*2.998e10)); //Update time
  Muon->ChangePosition(cpoint); //Update muon object position
  geom->GetCurrentTrack()->AddPoint(cpoint[0],cpoint[1],cpoint[2],Muon->GetTime()); //Add new point to the current track
  return;
}

//////////Photon Propagators//////////
//////////Photon Propagators//////////
//////////Photon Propagators//////////
//////////Photon Propagators//////////
//////////Photon Propagators//////////

void Tracker::Propagate_Photons(int n)
{
  int m = N_photons/n;

  for(int j=0;j<n;j++)
  {
    int i=m*j;
    InitializePhotonTrack(i);

    //Generate random step according to the probability of absorption of the photon
    //(the distance the photon goes through in the material before being absorbed)
    double absorption_step = generator->Generate_Photon_Step();
    double total_dist = 0; //Distance traveled by the photon in the material
    //int l =0;
    while(true)
    {
      //Find distance to the next boundary and set step
      geom->FindNextBoundary();

      if(CheckDensity()==1.023) //Photon is in scintillator
      {

        //If the absorption distance generated (minus the already travelled distance in this material) is shorter than the distance to the next boundary,
        //the photon is propagated until the point of absorption
        if((absorption_step - total_dist) < geom->GetStep())
        {
          geom->SetStep(absorption_step-total_dist);
          //Execute step (flag = kFALSE means it is an arbitrary step, not limited by geometrical reasons)
          geom->Step(kFALSE);
          Update_Photon_Track(i);

          N_absorbed++;
          break;
        } else {
          //Take step to next boundary but dont cross boundary (second flag is kFALSE)
          geom->Step(kTRUE, kFALSE);
          total_dist+=geom->GetStep();
          Update_Photon_Track(i);

          Photon_Scintillator_Reflection_Check(i);
        }

        if(DetectionCheck(cpoint,Photons[i]->GetEnergy()))
        {
          N_detected++;
          break;
        }

      }else
      {
        if(CheckDensity()==0) //Photon is in vacuum (or air gap)
        {
          //Take step to next boundary but dont cross boundary (second flag is kFALSE)
          geom->Step(kTRUE, kFALSE);
          Update_Photon_Track(i);

          if(!Photon_Vacuum_Reflection_Check(i)) //Photon is absorbed in the aluminium
          {
            N_lost++; //Photon lost to aluminium
            break;
          }

        }else //just in case the photon is inside the aluminium somehow
        {
          N_lost++;
          break;
        }
      }
    }

    delete Photons[i];
    //if(geom->IsOutside()){N_lost++;};

  }

  return;
}

void Tracker::InitializePhotonTrack(int i)
{

  //Add daughter track (photon track) to the main track and set it the current track
  geom->SetCurrentTrack(main_track->AddDaughter(i+1,Photons[i]->GetPDG(),Photons[i])); //id is set to i+1 because the main_track already has id=0
  //Add starting position to the track
  geom->GetCurrentTrack()->AddPoint(Photons[i]->GetStartingPosition()[0],Photons[i]->GetStartingPosition()[1],Photons[i]->GetStartingPosition()[2],Photons[i]->GetTime());
  //Setting both initial point and direction and finding the state
  geom->InitTrack(geom->GetCurrentTrack()->GetFirstPoint(), Photons[i]->GetDirection().data());
  //geom->SetCurrentPoint(Photons[i]->GetStartingPosition().data());
  return;
}

//Check if the photon in the scintillator is reflected in the boundary with the air gap
void Tracker::Photon_Scintillator_Reflection_Check(int i)
{
  vector<double> n = GetNormal(); //Get normal to the next boundary
  double theta = tools::Angle_Between_Vectors(Photons[i]->GetDirection(),n); //Compute incident angle

  if(CheckReflection(theta,1.58,1)) //Photon is reflected
  {
    vector<double> ndir = tools::Get_Reflected_Dir(Photons[i]->GetDirection(),n); //Get reflected direction
    geom->SetCurrentDirection(ndir.data()); //Set reflected direction the new direction
    Photons[i]->ChangeDirection(ndir); //Update photon object direction

  }else //Photon is transmited
  {
    vector<double> ndir = tools::Get_Refracted_Dir(Photons[i]->GetDirection(),n,theta,1.58,1); //Get refracted direction
    geom->FindNextBoundaryAndStep(); //Cross the boundary
    Update_Photon_Track(i);
    geom->SetCurrentDirection(ndir.data()); //Set refracted direction the new direction
    Photons[i]->ChangeDirection(ndir); //Update photon object direction

  }
  return;
}

//Check if the photon in the vacuum is reflected in the boundary with the scintillator or the aluminium foil
bool Tracker::Photon_Vacuum_Reflection_Check(int i)
{
  vector<double> n = GetNormal(); //Get normal to the next boundary
  double theta = tools::Angle_Between_Vectors(Photons[i]->GetDirection(),n); //Compute incident angle
  double r = sqrt(cpoint[0]*cpoint[0]+cpoint[1]*cpoint[1]); //Compute radial distance

  if(VacuumToPlastic(r)) //Photon is in the air gap near the boundary of the scintillator
  {
    if(CheckReflection(theta,1,1.58)) //Photon is reflected
    {
      vector<double> ndir = tools::Get_Reflected_Dir(Photons[i]->GetDirection(),n); //Get reflected direction
      geom->SetCurrentDirection(ndir.data()); //Set reflected direction the new direction
      Photons[i]->ChangeDirection(ndir); //Update photon object direction

    }else //Photon is transmited
    {
      vector<double> ndir = tools::Get_Refracted_Dir(Photons[i]->GetDirection(),n,theta,1,1.58); //Get refracted direction
      geom->FindNextBoundaryAndStep(); //Cross the boundary
      Update_Photon_Track(i);
      geom->SetCurrentDirection(ndir.data()); //Set refracted direction the new direction
      Photons[i]->ChangeDirection(ndir); //Update photon object direction
    }

  }else //Photon is near an aluminium foil boundary or near the top volume boundary
  {
    if(VacuumToAluminium(r)) //Photon is near an aluminium foil boundary
    {
      if(generator->Uniform(0,1) < .92) //Photon is reflected in the aluminium foil (with a probability of 92% for the relevant wavelengths)
      {
        vector<double> ndir = tools::Get_Reflected_Dir(Photons[i]->GetDirection(),n); //Get reflected direction
        geom->SetCurrentDirection(ndir.data()); //Set reflected direction the new direction
        Photons[i]->ChangeDirection(ndir); //Update photon object direction
        return true;
      } else{ //Photon is absorbed in the aluminium foil
        return false;
      }
    }else //Photon is near the top volume boundary (photon lost)
    {
      return false;
    }
  }

  return true;
}

void Tracker::Update_Photon_Track(int i){

  Photons[i]->IncreaseTime(GetRefractiveIndex()*geom->GetStep()/(2.998e10)); //Update time
  Photons[i]->ChangePosition(cpoint); //Update photon object position
  geom->GetCurrentTrack()->AddPoint(cpoint[0],cpoint[1],cpoint[2],Photons[i]->GetTime()); //Add new position to the track

}

//////////Draw Mode/////////
//////////Draw Mode/////////
//////////Draw Mode/////////
//////////Draw Mode/////////
//////////Draw Mode/////////

//Draw geometry and tracks
void Tracker::Draw(){
  //Create Canvas
  TCanvas *c1 = new TCanvas("c1","c1",1200,900);

  geom->GetTopVolume()->Draw(); //Draw geometry
  // /D: Track and first level descendents only are drawn
  // /*: Track and all descendents are drawn
  geom->DrawTracks("/*"); //Draw tracks
}

void Tracker::Heatmap()
{
  vector<double> phi;
  vector<double> z;
  int j=0;
  for(int i=0;i<N_photons;i++)
  {

    InitializePhotonTrack(i);
    geom->FindNextBoundary();
    geom->Step(kTRUE,kFALSE);
    if(cpoint[2]<0){break;};
    double r = sqrt(cpoint[0]*cpoint[0]+cpoint[1]*cpoint[1]);
    if(abs(r-Radius)<1e-6)
    {
      phi.push_back(tools::RadialTheta(cpoint));
      z.push_back(cpoint[2]-12.5);
      //cout<<phi[j]<<" "<<z[j]<<endl;
      //j++;
    }

  }
  auto H = new TH2D("h","",20,-M_PI,M_PI,20,0,1);
  for(int i=0;i<phi.size();i++)
  {
    H->Fill(phi[i],z[i]);
  }
  auto c = new TCanvas("c","c",1600,1000);
  H->Draw("colz");
  c->SaveAs("b.pdf");
  return;
}
