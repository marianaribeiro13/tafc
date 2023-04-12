#include "Tracker.h"

///////////////////////////// Constructor //////////////////////////////

Tracker::Tracker(double step, Generator* g,Geometry* G)
:  stepvalue(step), N_photons(0),N_absorbed(0),N_detected(0),N_lost(0)
{
  generator = g;
  Geo = G;
  //Build the Telescope Geometry
  nav = Geo->GetGeoManager()->GetCurrentNavigator();
  if (!nav) nav = Geo->GetGeoManager()->AddNavigator();


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
    double n_density=3.33e29;
    double I = 1.03660828e-17;
    return ((qe*qe*qe*qe*n_density*Z*Z*(log((2*me*c*c*x[0]*x[0])/(I*(1-x[0]*x[0])))-x[0]*x[0]))/(4*M_PI*eps0*eps0*me*c*c*x[0]*x[0]))/(1.602e-13);
  };
  BetheBloch = new TF1("f",f); //Create BetheBloch TF1

  //Generate Cosmic Muon at incident scintillator plane
  Muon = generator->Generate_CosmicMuon(generator->Generate_Position(Geo->GetDistance(),Geo->GetHeight(),Geo->GetRadius()));
  MainTrackID = Geo->GetGeoManager()->AddTrack(0,Muon->GetPDG(),Muon);
  main_track = Geo->GetGeoManager()->GetTrack(MainTrackID);
  main_track->AddPoint(Muon->GetStartingPosition()[0],Muon->GetStartingPosition()[1],Muon->GetStartingPosition()[2],0);
  Geo->GetGeoManager()->SetCurrentTrack(main_track);
  //Add muon track to the Navigator and assign it to main_track atribute
  nav->InitTrack(Muon->GetStartingPosition().data(), Muon->GetDirection().data());

  //Get pointer to current position (Whenever we update the position cpoint is updated,
  //since it is equal to the nav pointer to the current position)
  cpoint = nav->GetCurrentPoint();


  Photons_flag = false; // This Flag Checks if the photons have been propagated or not
  //MainTrackID = 0;      //The ID of the main muon track
  DoubleCross = false;  //This Flag Checks if the muon has crossed both detectors
}

///////////////////////////////// Destructor ///////////////////////////////////

Tracker::~Tracker()
{
  if(!Photons_flag)
  {
    Photons.clear();
  }

  delete Muon;
  main_track->ResetTrack();

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
  return nav->IsSameLocation(aux[0],aux[1],aux[2]);
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

//Check if light is reflected or transmitted and get new photon direction
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
  vector<double> d = {nav->GetCurrentDirection()[0],nav->GetCurrentDirection()[1],nav->GetCurrentDirection()[2]};



  return Geo->GetNormal(r,h,d,cpoint);
}

//Check if the particle is in vacuum near the scintillator boundary
// bool Tracker::VacuumToPlastic(double r)
// {
//   //If the particle is in vacuum (in the air gap) and near the lateral scintillator boundary
//   //or near one of the flat scintillator boundaries
//   if(CheckDensity()==0 && ((abs(r-Radius) < 1e-6) || (abs(abs(cpoint[2])-(0.5*Distance+Height))<1e-6)
//   || (abs(abs(cpoint[2])-(0.5*Distance))<1e-6))){return true;};
//
//   return false;
// }

//Check if the particle is in vacuum near aluminium foil
// bool Tracker::VacuumToAluminium(double r)
// {
//   double h = abs(abs(cpoint[2])-0.5*(Distance+Height));
//
//   if(CheckDensity()==0 && ((abs(r-innerradius) < 1e-6) || (abs(r-outerradius) <1e-6)
//   ||abs(h-(Height/2+Airgap))<1e-6 )){return true;};
//
//   return false;
// }

//Check if photon is in the detector region and if it is detected according to the detector efficiency
bool Tracker::DetectionCheck(const double *cpoint,double e)
{
  if(Geo->Check_Symmetric_Detector(cpoint) && generator->Uniform(0,1)<generator->GetDetector_Efficiency()->Eval(e))
  {
    return true;
  }
  return false;
}

double Tracker::CheckDensity()
{
    return nav->GetCurrentNode()->GetVolume()->GetMedium()->GetMaterial()->GetDensity();
}

double Tracker::GetRefractiveIndex()
{
  if(CheckDensity() == 1.023)
  {
    return 1.58;
  }
  if(CheckDensity()==0)
  {
    return 1;
  }else
  {
    return 1.58;
  }
  return 0;
}

//////////Simulation Control//////////
//////////Simulation Control//////////
//////////Simulation Control//////////
//////////Simulation Control//////////
//////////Simulation Control//////////

//Deletes the muon and all the photons, resets the counters and flags and clears all tracks from the geometry
void Tracker::Reset()
{
  if(!Photons_flag)
  {
    for(int i =0;i<Photons.size();i++)
    {
      delete Photons[i];
    }
  }
  main_track->ResetTrack();
  Geo->GetGeoManager()->SetCurrentTrack(main_track);
  Photons.clear();
  delete Muon;
  Photons_flag = false;
  DoubleCross = false;
  N_photons = 0;
  N_absorbed = 0;
  N_detected = 0;
  N_lost = 0;
  return;
}

//Resets the tracker then andomly generates a new muon on the scintillator plane and initializes the main track
void Tracker::Generate_New_Muon()
{
  Reset();
  Muon = generator->Generate_CosmicMuon(generator->Generate_Position(Geo->GetDistance(),Geo->GetHeight(),Geo->GetRadius()));

  main_track->AddPoint(Muon->GetStartingPosition()[0],Muon->GetStartingPosition()[1],Muon->GetStartingPosition()[2],0);
  nav->InitTrack(main_track->GetFirstPoint(), Muon->GetDirection().data());
  cpoint = nav->GetCurrentPoint();

  return;
}

//Resets the tracker then inserts a manually generated muon and initializes the main track
void Tracker::Insert_New_Muon(Particle* M)
{
  Reset();
  Muon = M;

  main_track->AddPoint(Muon->GetStartingPosition()[0],Muon->GetStartingPosition()[1],Muon->GetStartingPosition()[2],0);
  nav->InitTrack(main_track->GetFirstPoint(), Muon->GetDirection().data());
  cpoint = nav->GetCurrentPoint();


  return;
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

  while(!nav->IsOutside()) //While the particle is inside the defined top volume
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
  nav->FindNextBoundaryAndStep(); //Make step to the next boundary and cross it (updates current position of the navigator)
  Muon->IncreaseTime(nav->GetStep()/(Muon->GetVelocity()*2.998e10)); //Update time
  Muon->ChangePosition(cpoint); //Update muon object position
  main_track->AddPoint(cpoint[0],cpoint[1],cpoint[2],Muon->GetTime()); //Add new point to the current track
  return;
}

//Make scintillator steps until the particle leaves the scintillator
void Tracker::Muon_Scintillator_Step()
{
  //Set next step to the defined arbitrary step value (defined by the constructer)
  nav->SetStep(stepvalue);
  while(CheckSameLocation()) //While the next step is still inside the same volume
  {
    //Execute step defined by SetStep (flag = kFALSE means it is an arbitrary step, not limited by geometrical reasons)
    nav->Step(kFALSE); //Updates the current position of the navigator
    Muon->IncreaseTime(stepvalue/(Muon->GetVelocity()*2.998e10)); //Update time
    Muon->ChangePosition(cpoint); //Update muon object position
    main_track->AddPoint(cpoint[0],cpoint[1],cpoint[2],Muon->GetTime()); //Add new point to the current track

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

  nav->FindNextBoundaryAndStep(stepvalue); //Make step to the next boundary and cross it (step is less than the defined step)
  Muon->IncreaseTime(nav->GetStep()/(Muon->GetVelocity()*2.998e10)); //Update time
  Muon->ChangePosition(cpoint); //Update muon object position
  main_track->AddPoint(cpoint[0],cpoint[1],cpoint[2],Muon->GetTime()); //Add new point to the current track

  //Get number of photons from Poisson distribution (usually approximatly 1 photon per 100 eV for common scintillators)
  int n = generator->Generate_Photon_Number(10000*Update_Energy(nav->GetStep()));
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
  nav->FindNextBoundaryAndStep(); //Make step to the next boundary and cross it (updates current position of the navigator)
  Muon->IncreaseTime(nav->GetStep()/(Muon->GetVelocity()*2.998e10)); //Update time
  Muon->ChangePosition(cpoint); //Update muon object position
  main_track->AddPoint(cpoint[0],cpoint[1],cpoint[2],Muon->GetTime()); //Add new point to the current track
  return;
}

//////////Photon Propagators//////////
//////////Photon Propagators//////////
//////////Photon Propagators//////////
//////////Photon Propagators//////////
//////////Photon Propagators//////////

void Tracker::Propagate_Photons()
{

  Photons_flag = true;
  for(int i=0;i<N_photons;i++)
  {

    InitializePhotonTrack(i);

    //Generate random step according to the probability of absorption of the photon
    //(the distance the photon goes through in the material before being absorbed)
    double absorption_step = generator->Generate_Photon_Step();
    double total_dist = 0; //Distance traveled by the photon in the material

    while(true)
    {
      //Find distance to the next boundary and set step
      nav->FindNextBoundary();
      nav->Step(kTRUE,kFALSE);

      if(CheckDensity()==1.023)
      {

        if((total_dist+=nav->GetStep()) >=absorption_step)
        {
          N_absorbed++;
          break;
        }

        Photon_Scintillator_Reflection_Check(i);

        if(DetectionCheck(cpoint,Photons[i]->GetEnergy()))
        {
          N_detected++;
          break;
        }
      }else
      {
        if(CheckDensity()==0)
        {

          if(!Photon_Vacuum_Reflection_Check(i))
          {
            N_lost++; //Photon lost to aluminium
            break;
          }
        }else
        {
          N_lost++;
          break;
        }
      }
      Update_Photon_Track(i);
      // Photons[i]->IncreaseTime(GetRefractiveIndex()*nav->GetStep()/(2.998e10));
      // Photons[i]->ChangePosition(cpoint);
      // nav->GetCurrentTrack()->AddPoint(cpoint[0],cpoint[1],cpoint[2],Photons[i]->GetTime());
    }
    delete Photons[i];
  }
  return;
}

void Tracker::InitializePhotonTrack(int i)
{
  //Add daughter track (photon track) to the main track and set it the current track
  current_track = main_track->AddDaughter(i+1,Photons[i]->GetPDG(),Photons[i]);
  Geo->GetGeoManager()->SetCurrentTrack(current_track); //id is set to i+1 because the main_track already has id=0
  //Add starting position to the track
  current_track->AddPoint(Photons[i]->GetStartingPosition()[0],Photons[i]->GetStartingPosition()[1],Photons[i]->GetStartingPosition()[2],Photons[i]->GetTime());
  //Setting both initial point and direction and finding the state
  nav->InitTrack(Geo->GetGeoManager()->GetCurrentTrack()->GetFirstPoint(), Photons[i]->GetDirection().data());
  //geom->SetCurrentPoint(Photons[i]->GetStartingPosition().data());
  return;
}

void Tracker::Photon_Scintillator_Reflection_Check(int i)
{
  vector<double> n = GetNormal();
  double theta = tools::Angle_Between_Vectors(Photons[i]->GetDirection(),n);
  if(CheckReflection(theta,1.58,1) )
  {
    vector<double> ndir = tools::Get_Reflected_Dir(Photons[i]->GetDirection(),n);
    nav->SetCurrentDirection(ndir.data());
    Photons[i]->ChangeDirection(ndir);

  }else
  {
    vector<double> ndir = tools::Get_Refracted_Dir(Photons[i]->GetDirection(),n,theta,1.58,1);
    nav->FindNextBoundaryAndStep();
    nav->SetCurrentDirection(ndir.data());
    Photons[i]->ChangeDirection(ndir);

  }
  return;
}

bool Tracker::Photon_Vacuum_Reflection_Check(int i)
{

  vector<double> n = GetNormal();
  double theta = tools::Angle_Between_Vectors(Photons[i]->GetDirection(),n);
  double r = sqrt(cpoint[0]*cpoint[0]+cpoint[1]*cpoint[1]);
  double h = abs(abs(cpoint[2])-0.5*(Geo->GetDistance()+Geo->GetHeight()));
  if(Geo->VacuumToPlastic(r,h))
  {
    if(CheckReflection(theta,1,1.58) )
    {
      vector<double> ndir = tools::Get_Reflected_Dir(Photons[i]->GetDirection(),n);
      nav->SetCurrentDirection(ndir.data());
      Photons[i]->ChangeDirection(ndir);
    }else
    {
      vector<double> ndir = tools::Get_Refracted_Dir(Photons[i]->GetDirection(),n,theta,1,1.58);
      nav->SetCurrentDirection(ndir.data());
      Photons[i]->ChangeDirection(ndir);
      nav->FindNextBoundaryAndStep();
    }
  }else
  {
    if(Geo->VacuumToAluminium(r,h))
    {
      if(generator->Uniform(0,1) < .92)
      {
        vector<double> ndir = tools::Get_Reflected_Dir(Photons[i]->GetDirection(),n);
        nav->SetCurrentDirection(ndir.data());
        Photons[i]->ChangeDirection(ndir);
        return true;
      }
    }else
    {
      return false;
    }
  }
  return true;
}

void Tracker::Update_Photon_Track(int i){

  Photons[i]->IncreaseTime(GetRefractiveIndex()*nav->GetStep()/(2.998e10)); //Update time
  Photons[i]->ChangePosition(cpoint); //Update photon object position
  current_track->AddPoint(cpoint[0],cpoint[1],cpoint[2],Photons[i]->GetTime()); //Add new position to the track
  return;
}

//////////Draw Mode/////////
//////////Draw Mode/////////
//////////Draw Mode/////////
//////////Draw Mode/////////
//////////Draw Mode/////////

//Propagate and draw n photon tracks
void Tracker::Draw(int n)
{

  int m = N_photons/n; //Only propagate the mth photon, in total propagating n photons,
  for(int j=0;j<n;j++)
  {
    int i=m*j;
    InitializePhotonTrack(i);

    double absorption_step = generator->Generate_Photon_Step();
    double total_dist = 0;

    while(true) //This part is the same as the photon propagator
    {

      nav->FindNextBoundary();
      nav->Step(kTRUE,kFALSE);

      if(CheckDensity()==1.023)
      {

        if((total_dist+=nav->GetStep()) >=absorption_step)
        {
          N_absorbed++;
          break;
        }

        Photon_Scintillator_Reflection_Check(i);

        if(DetectionCheck(cpoint,Photons[i]->GetEnergy()))
        {
          N_detected++;
          break;
        }
      }else
      {
        if(CheckDensity()==0)
        {

          if(!Photon_Vacuum_Reflection_Check(i))
          {
            N_lost++; //Photon lost to aluminium
            break;
          }
        }else
        {
          N_lost++;
          break;
        }
      }

      Update_Photon_Track(i);
    }
  }

  TCanvas *c1 = new TCanvas("c1","c1",1200,900);
  Geo->GetGeoManager()->GetTopVolume()->Draw(); //Draw geometry
  // /D: Track and first level descendents only are drawn
  // /*: Track and all descendents are drawn
  Geo->GetGeoManager()->DrawTracks("/*"); //Draw tracks
}

//Fills the given Histogram with the detection points of all photons in the first scintillator
void Tracker::Fill_Heatmap(TH2D* H)
{

  for(int i=0;i<N_photons;i++)
  {

    InitializePhotonTrack(i);
    nav->FindNextBoundary();
    nav->Step(kTRUE,kFALSE);
    if(cpoint[2]>0)
    {
      double r = sqrt(cpoint[0]*cpoint[0]+cpoint[1]*cpoint[1]);
      if(abs(r-Geo->GetRadius())<1e-6)
      {
        H->Fill(tools::RadialTheta(cpoint),cpoint[2]-12.5);
      }
    }



  }
  return;
}
