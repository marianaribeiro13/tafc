#include "Tracker.h"
#include <mutex>

//#define DRAWMODE

#ifdef DRAWMODE
std::mutex Mutex;
#endif 
///////////////////////////// Constructer //////////////////////////////

Tracker::Tracker(TGeoManager* GeoM, Generator* g, Particle* part)
: N_photons(0), N_absorbed(0),N_detected(0), N_photons1(0), N_lost(0), DoubleCross(false)
{
  stepvalue = param.step;

  geom = GeoM; //Assign pointer to TGeoManager object
  generator = g; //Assign pointer to Generator object
  Muon = part; //Assign pointer to Particle object
  main_track = nullptr; //Just in case draw mode is not defined (not sure if this should be done)

  //The functions for adding navigators have a lock inside (they are thread safe)
  //Get the current navigator for this thread (if it was already added 
  //- happens if we create more than one Tracker object in the same thread)
  nav = geom->GetCurrentNavigator(); 
  if (!nav) nav = geom->AddNavigator();

#ifdef DRAWMODE
  Mutex.lock();
  //Add muon track to the Navigator and assign it to main_track atribute
  int track_index = geom->AddTrack(0,Muon->GetPDG(),Muon); //The first argument is not important
  Mutex.unlock();

  main_track = geom->GetTrack(track_index); //Get muon track
  //Add starting point to the track (initial position of the muon) and set main_track as the current track
  main_track->AddPoint(Muon->GetStartingPosition()[0],Muon->GetStartingPosition()[1],Muon->GetStartingPosition()[2],0);
  geom->SetCurrentTrack(track_index);
 
#endif
  
  //Get pointer to current position (Whenever we update the position cpoint is updated, 
  //since it is equal to the nav pointer to the current position)
  cpoint = nav->GetCurrentPoint(); 
  //Get pointer to current direction (Whenever we update the direction cdir is updated, 
  //since it is equal to the nav pointer to the current direction)
  cdir = nav->GetCurrentDirection();

  Photons_flag = false; // This Flag Checks if the photons have been propagated or not
  //MainTrackID = 0;      //The ID of the main muon track
  DoubleCross = false;  //This Flag Checks if the muon has crossed both scintillators
}

///////////////////////////////// Destructor ///////////////////////////////////

Tracker::~Tracker()
{
  if(Muon){
    delete Muon;
  }
  //if(!Photons_flag)
  //{
  for(int i =0;i<Photons.size();i++){
    if(Photons[i]){
      delete Photons[i];
    }
  }
  
  //if(cpoint){
    //delete cpoint;
  //}

  //if(cdir){
    //delete cdir;
  //}
  //}
}

//Update energy and momentum of the particle after energy loss by BetheBloch equation
double Tracker::Update_Energy(double step)
{
  double dE = BetheBloch(Muon->GetVelocity()) * (step/100); //step converted from cm to m
  Muon->ChangeEnergy(Muon->GetEnergy()-dE); //Update energy
  Muon->ChangeMomentum(Muon->CalculateMomentum(Muon->GetEnergy())); //Update momentum
  return dE;

}

//Checks whether the next position (if we make the defined step value in the current particle direction)
//is in the same node/volume of the current position
bool Tracker::CheckSameLocation() //Parallel
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
  //vector<double> d = {cdir[0], cdir[1], cdir[2]};
  vector<double> d = {nav->GetCurrentDirection()[0],nav->GetCurrentDirection()[1],nav->GetCurrentDirection()[2]};

  if(abs(r-param.Radius) <1e-6 || abs(r-param.innerradius) <1e-6 || abs(r-param.outerradius) <1e-6)
  {

    vertical_reflection=true;
    aux[0] = cpoint[0]/r;
    aux[1] = cpoint[1]/r;
    aux[2] = 0;
  }
  if((abs(h-(0.5*param.Distance+param.Height))<1e-6) || (abs(h-0.5*param.Distance)<1e-6) || (abs(h-(param.Airgap+(0.5*param.Distance)+param.Height))<1e-6)  || (abs(h-(-param.Airgap+(0.5*param.Distance)))<1e-6) || (abs(h-(param.Thickness+param.Airgap+(0.5*param.Distance)+param.Height))<1e-6) || (abs(h-(-param.Thickness-param.Airgap+(0.5*param.Distance)))<1e-6) )
  {
    horizontal_reflection =true;
    aux[0] = 0;
    aux[1] = 0;
    aux[2] = 1;
  }
  if(vertical_reflection && horizontal_reflection)// Corner reflection, photon goes back
  {
    return d;
  }

  if(tools::Angle_Between_Vectors(aux,d)>M_PI/2)
  {
    for(int i=0;i<3;i++){aux[i] = -aux[i];};
  }

  return aux;
}

//Check if the particle is in vacuum near the scintillator boundary
bool Tracker::VacuumToPlastic(double r)
{
  double h = abs(abs(cpoint[2])-0.5*(param.Distance+param.Height));
  //If the particle is in vacuum (in the air gap) and near the lateral scintillator boundary
  //or near one of the flat scintillator boundaries
  if(((abs(r-param.Radius) < 1e-6) || abs(h-(param.Height/2))<1e-6)){return true;};

  return false;

  
}

//Check if the particle is in vacuum near aluminium foil
bool Tracker::VacuumToAluminium(double r)
{
  double h = abs(abs(cpoint[2])-0.5*(param.Distance+param.Height));
  //If the particle is in vacuum (either in the air gap or outside the telescope) and near any aluminium foil boundary
  if(((abs(r-param.innerradius) < 1e-6) || (abs(r-param.outerradius) <1e-6) || abs(h-(param.Height/2+param.Airgap))<1e-6 )){return true;};

  return false;
}

//Check if photon is in the detector region and if it is detected according to the detector efficiency
bool Tracker::DetectionCheck(double e)
{

  if(Check_Symmetric_Detector() && generator->Uniform(0,1)<generator->GetDetector_Efficiency()->Eval(e))
  {

    return true;
  }

  return false;
}

//Check if photon is in the detector region and if it is detected according to the detector efficiency
bool Tracker::Is_Photon_Detected(double E)
{

  if(Is_Detector_Region() && generator->Uniform(0,1)<generator->GetDetector_Efficiency()->Eval(E))
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

  //Setting both initial point and direction and finding the state
  nav->InitTrack(Muon->GetStartingPosition().data(), Muon->GetDirection().data());

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
      if(scintillator_cross == 1){
        N_photons1 = N_photons;
      }
    }
    if(CheckDensity()==2.7) //Particle is in the aluminium foil
    {
      Muon_Aluminium_Step();
    }
  }

  if(scintillator_cross < 2){ //The particle crosses both scintillators
#ifdef DRAWMODE
    //Mutex.lock();
    //geom->ClearTracks();
    main_track->ResetTrack();
    //Mutex.unlock();
#endif
  }else { 
    DoubleCross=true;
  }

  return;
}

//Make vacuum step
void Tracker::Muon_Vacuum_Step()
{
  nav->FindNextBoundaryAndStep(); //Make step to the next boundary and cross it (updates current position of the navigator)
  Muon->ChangePosition(cpoint); //Update muon object position
#ifdef DRAWMODE
  //Mutex.lock();
  Muon->IncreaseTime(nav->GetStep()/(Muon->GetVelocity()*2.998e10)); //Update time
  main_track->AddPoint(cpoint[0],cpoint[1],cpoint[2],Muon->GetTime()); //Add new point to the current track
  //Mutex.unlock();
#endif
  return;
}

//Make scintillator steps until the particle leaves the scintillator
void Tracker::Muon_Scintillator_Step()
{
  //Set next step to the defined arbitrary step value (defined by the constructer)
  nav->SetStep(stepvalue);
  while(CheckSameLocation())
  {
    //Execute step defined by SetStep (flag = kFALSE means it is an arbitrary step, not limited by geometrical reasons)
    nav->Step(kFALSE); //Updates the current position of the navigator
    Muon->ChangePosition(cpoint); //Update muon object position
#ifdef DRAWMODE
    //Mutex.lock();
    Muon->IncreaseTime(stepvalue/(Muon->GetVelocity()*2.998e10)); //Update time
    main_track->AddPoint(cpoint[0],cpoint[1],cpoint[2],Muon->GetTime()); //Add new point to the current track
    //Mutex.unlock();
#endif
    //Get number of photons from Poisson distribution (usually approximatly 1 photon per 100 eV for common scintillators)
    int n = generator->Generate_Photon_Number(10000*Update_Energy(stepvalue)); //Energy is converted to MeV
    N_photons+=n; //Update total number of emited photons

    for(int i=0;i<n;i++)
    {
      //Add generated photon according to Scintillator spectrum to the vector of photons (to propagate later)
      Photons.emplace_back(generator->Generate_Photon(Muon->GetPosition())); 
      Photons[i]->IncreaseTime(Muon->GetTime()); //Assign current time to the photon
    }
  }

  nav->FindNextBoundaryAndStep(stepvalue); //Make step to the next boundary and cross it (step is less than the defined step)
  Muon->ChangePosition(cpoint); //Update muon object position
#ifdef DRAWMODE
  //Mutex.lock();
  Muon->IncreaseTime(nav->GetStep()/(Muon->GetVelocity()*2.998e10)); //Update time
  main_track->AddPoint(cpoint[0],cpoint[1],cpoint[2],Muon->GetTime()); //Add new point to the current track
  //Mutex.unlock();
#endif
  //Get number of photons from Poisson distribution (usually approximatly 1 photon per 100 eV for common scintillators)
  int n = generator->Generate_Photon_Number(10000*Update_Energy(nav->GetStep()));
  N_photons+=n; //Update total number of emited photons

  for(int i=0;i<n;i++)
  {
    //Add generated photon according to Scintillator spectrum to the vector of photons (to propagate later)
    Photons.emplace_back(generator->Generate_Photon(Muon->GetPosition()));
    Photons[i]->IncreaseTime(Muon->GetTime()); //Assign current time to the photon
  }

  return;
}

//Make aluminium step
void Tracker::Muon_Aluminium_Step()
{
  nav->FindNextBoundaryAndStep(); //Make step to the next boundary and cross it (updates current position of the navigator)
  Muon->ChangePosition(cpoint); //Update muon object position
#ifdef DRAWMODE
  //Mutex.lock();
  Muon->IncreaseTime(nav->GetStep()/(Muon->GetVelocity()*2.998e10)); //Update time
  main_track->AddPoint(cpoint[0],cpoint[1],cpoint[2],Muon->GetTime()); //Add new point to the current track
  //Mutex.unlock();
#endif
  return;
}

//////////Photon Propagators//////////
//////////Photon Propagators//////////
//////////Photon Propagators//////////
//////////Photon Propagators//////////
//////////Photon Propagators//////////

void Tracker::Propagate_Photons(int iphoton, int fphoton)
{
  Photons_flag = true;

  //int m = N_photons/n;
  for(int i=iphoton; i<fphoton;i++)
  {
    //int i=m*j;
    nav->InitTrack(Photons[i]->GetStartingPosition().data(), Photons[i]->GetDirection().data());
#ifdef DRAWMODE
    //Mutex.lock();
    InitializePhotonTrack(i);
    //Mutex.unlock();
#endif
    //Generate random step according to the probability of absorption of the photon 
    //(the distance the photon goes through in the material before being absorbed)
    double absorption_step = generator->Generate_Photon_Step();
    double total_dist = 0; //Distance traveled by the photon in the material
    //int l =0;
    while(true)
    {
      //Find distance to the next boundary and set step
      nav->FindNextBoundary();

      if(CheckDensity()==1.023) //Photon is in scintillator
      {
        
        //If the absorption distance generated (minus the already travelled distance in this material) is shorter than the distance to the next boundary, 
        //the photon is propagated until the point of absorption
        if((absorption_step - total_dist) < nav->GetStep())
        {
          nav->SetStep(absorption_step-total_dist);
          //Execute step (flag = kFALSE means it is an arbitrary step, not limited by geometrical reasons)
          nav->Step(kFALSE);
          Update_Photon(i);
          N_absorbed++;
          break;
        } else {
          //Take step to next boundary but dont cross boundary (second flag is kFALSE)
          nav->Step(kTRUE, kFALSE);
          total_dist += nav->GetStep();
          Update_Photon(i);

          Photon_Scintillator_Reflection_Check(i);
        }

        if(DetectionCheck(Photons[i]->GetEnergy())) //Check if photon is detected
        {
          N_detected++;
          break;
        }

      }else
      {
        if(CheckDensity()==0) //Photon is in vacuum (or air gap)
        {
          //Take step to next boundary but dont cross boundary (second flag is kFALSE)
          nav->Step(kTRUE, kFALSE);
          Update_Photon(i);

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

    //delete Photons[i];
    //std::cout << Photons.size() << '\n';

  }

  return;
}

void Tracker::InitializePhotonTrack(int i)
{
#ifdef DRAWMODE
  Mutex.lock();
  //Add daughter track (photon track) to the main track and set it the current track
  geom->SetCurrentTrack(main_track->AddDaughter(0,Photons[i]->GetPDG(),Photons[i])); //the first argument is not important
  Mutex.unlock();
  //Add starting position to the track
  geom->GetCurrentTrack()->AddPoint(Photons[i]->GetStartingPosition()[0],Photons[i]->GetStartingPosition()[1],Photons[i]->GetStartingPosition()[2],Photons[i]->GetTime());
#endif
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
    nav->SetCurrentDirection(ndir.data()); //Set reflected direction the new direction
    Photons[i]->ChangeDirection(ndir); //Update photon object direction

  }else //Photon is transmited
  {
    vector<double> ndir = tools::Get_Refracted_Dir(Photons[i]->GetDirection(),n,theta,1.58,1); //Get refracted direction
    nav->FindNextBoundaryAndStep(); //Cross the boundary
    Update_Photon(i);
    nav->SetCurrentDirection(ndir.data()); //Set refracted direction the new direction
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
      nav->SetCurrentDirection(ndir.data()); //Set reflected direction the new direction
      Photons[i]->ChangeDirection(ndir); //Update photon object direction
    
    }else //Photon is transmited
    {
      vector<double> ndir = tools::Get_Refracted_Dir(Photons[i]->GetDirection(),n,theta,1,1.58); //Get refracted direction
      nav->FindNextBoundaryAndStep(); //Cross the boundary
      Update_Photon(i);
      nav->SetCurrentDirection(ndir.data()); //Set refracted direction the new direction
      Photons[i]->ChangeDirection(ndir); //Update photon object direction
    }

  }else //Photon is near an aluminium foil boundary or near the top volume boundary
  {
    if(VacuumToAluminium(r)) //Photon is near an aluminium foil boundary
    {
      if(generator->Uniform(0,1) < .92) //Photon is reflected in the aluminium foil (with a probability of 92% for the relevant wavelengths)
      {
        vector<double> ndir = tools::Get_Reflected_Dir(Photons[i]->GetDirection(),n); //Get reflected direction
        nav->SetCurrentDirection(ndir.data()); //Set reflected direction the new direction
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

void Tracker::Update_Photon(int i){

  Photons[i]->ChangePosition(cpoint); //Update photon object position
#ifdef DRAWMODE
  //Mutex.lock();
  Photons[i]->IncreaseTime(GetRefractiveIndex()*nav->GetStep()/(2.998e10)); //Update time
  geom->GetCurrentTrack()->AddPoint(cpoint[0],cpoint[1],cpoint[2],Photons[i]->GetTime()); //Add new position to the track
  //Mutex.unlock();
#endif
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

bool Tracker::Check_Symmetric_Detector()
{
  double h = abs(cpoint[2]);

  if(h>(0.5*(param.Height+param.Distance+param.SIPM_size)) || h<(0.5*(param.Height+param.Distance-param.SIPM_size))){return false;};

  double r = sqrt(cpoint[0]*cpoint[0]+cpoint[1]*cpoint[1]);
  if(abs(r-param.Radius)>1e-6){return false;};

  double theta = atan(cpoint[1]/cpoint[0])/param.SIPM_angle;
  double delta = theta - round(theta);
  if(abs(delta) > param.SIPM_alpha){return false;};

  return true;
}

bool Tracker::Is_Detector_Region()
{
  //Check if the photon is in the right z range
  double h = abs(cpoint[2]);
  if(h>(0.5*(param.Height+param.Distance+param.SIPM_size)) || h<(0.5*(param.Height+param.Distance-param.SIPM_size))){return false;};

  //Check if the photon is in the scintillator boundary
  double r = sqrt(cpoint[0]*cpoint[0]+cpoint[1]*cpoint[1]);
  if(abs(r-param.Radius)>1e-6){return false;};

  //Current phi angle of the photon
  double phi = tools::PhiAngle(cpoint);
  
  //Check if the photon is in the phi angular range of any SIPM
  for(double SIPM_phi : param.SIPM_angles){
    if(abs(phi - SIPM_phi) < param.SIPM_phi_range/2){
      return true;
    }
  }

  return false;
}

////////////////////////////////////// Bethe Bloch function to calculate energy loss in material //////////////////////////////////

double Tracker::BetheBloch(double v){

    int Z=-1;
    double c = 299792458;
    double mp = 1.672621637e-27;
    double me = 9.1093821499999992e-31;
    double qe = 1.602176487e-19;
    double na = 6.02214179e23;
    double eps0 = 8.854187817e-12;
    double n_density=3.33e29; // per m^3
    double I = 1.03660828e-17;

    return ((qe*qe*qe*qe*n_density*Z*Z*(log((2*me*c*c*v*v)/(I*(1-v*v)))-v*v))/(4*M_PI*eps0*eps0*me*c*c*v*v))/(1.602e-13);
}

std::vector<std::pair<double,double>> Tracker::PropagatePhotons_To_FirstBoundary(int iphoton, int fphoton)
{
  std::vector<std::pair<double,double>> points;
  
  for(int i=iphoton; i<fphoton; i++){

    nav->InitTrack(Photons[i]->GetStartingPosition().data(), Photons[i]->GetDirection().data());
    nav->FindNextBoundary();
    nav->Step(kTRUE,kFALSE);
    double r = sqrt(cpoint[0]*cpoint[0]+cpoint[1]*cpoint[1]);
    
    if(abs(r-param.Radius)<1e-6){
      points.emplace_back(std::make_pair(tools::PhiAngle(cpoint), cpoint[2] - param.Distance/2));
    }
  }

  return points;
}


//Deletes the muon and all the photons, resets the counters and flags and clears all tracks from the geometry
/*void Tracker::Reset(){
  
#ifdef DRAWMODE
  //geom->ClearTracks();
  main_track->ResetTrack();
#endif
}*/