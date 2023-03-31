#include "Tracker.h"

Tracker::Tracker(double radius, double height, double distance, double airgap, double althickness, double step, Generator* g, int n_SIPMS)
: Geometry(), stepvalue(step), N_photons(0),N_absorbed(0),N_detected(0),N_lost(0),DoubleCross(false)
{
  generator = g;
  Build_MuonTelescope(radius, height, distance, airgap, althickness,n_SIPMS);
  Muon = generator->Generate_CosmicMuon(generator->Generate_Position(Distance,Height,Radius));

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
  BetheBloch = new TF1("f",f);
  geom->AddTrack(0,Muon->GetPDG(),Muon);
  geom->GetTrack(0)->AddPoint(Muon->GetStartingPosition()[0],Muon->GetStartingPosition()[1],Muon->GetStartingPosition()[2],0);
  geom->SetCurrentTrack(0);
  geom->InitTrack(geom->GetCurrentTrack()->GetFirstPoint(), Muon->GetDirection().data());
  cpoint = geom->GetCurrentPoint();
}

Tracker::~Tracker(){}

double Tracker::CheckDensity()
{
    return geom->GetCurrentNode()->GetVolume()->GetMedium()->GetMaterial()->GetDensity();
}

double Tracker::Update_Energy(double step)
{
  double dE = BetheBloch->Eval(Muon->GetVelocity()) * (step/100);
  Muon->ChangeEnergy(Muon->GetEnergy()-dE);
  Muon->ChangeMomentum(Muon->CalculateMomentum(Muon->GetEnergy()));
  return dE;

}

bool Tracker::CheckSameLocation()
{
  double aux[3];
  for(int i=0;i<3;i++)
  {
    aux[i] = cpoint[i] + stepvalue*Muon->GetDirection()[i];
  }
  return geom->IsSameLocation(aux[0],aux[1],aux[2]);
}

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

bool Tracker::CheckReflection(double thetai,double n1,double n2){


    if(n1>n2 && thetai > asin(n2/n1))
    {
      return true;
    }

    double Reff = FresnelLaw(thetai, n1, n2);

    if(generator->Uniform(0,1) < Reff) {

        return true;

    } else {

        return false;
    }
}

vector<double> Tracker::GetNormal()
{
  double r = sqrt(cpoint[0]*cpoint[0]+cpoint[1]*cpoint[1]);
  double h = abs(cpoint[2]);
  vector<double> aux(3);
  bool horizontal_reflection =false,vertical_reflection=false;
  vector<double> d = {geom->GetCurrentDirection()[0],geom->GetCurrentDirection()[1],geom->GetCurrentDirection()[2]};

  if(abs(r-Radius) <1e-6 || abs(r-innerradius) <1e-6 || abs(r-outerradius) <1e-6)
  {

    vertical_reflection=true;
    aux[0] = cpoint[0]/r;
    aux[1] = cpoint[1]/r;
    aux[2] = 0;
  }
  if((abs(h-(0.5*Distance+Height))<1e-6) || (abs(h-0.5*Distance)<1e-6) || (abs(h-(Airgap+(0.5*Distance+Height)))<1e-6)  || (abs(h-(-Airgap+(0.5*Distance)))<1e-6) || (abs(h-(Thickness+Airgap+(0.5*Distance+Height)))<1e-6) || (abs(h-(-Thickness-Airgap+(0.5*Distance)))<1e-6)  )
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

bool Tracker::ReflectionHandler(int i)
{
  double n1=1,n2=1;
  double r = sqrt(cpoint[0]*cpoint[0]+cpoint[1]*cpoint[1]);
  vector<double> n = GetNormal();

  if(CheckDensity() == 1.023)
  {
    n1 = 1.58;
    n2 = 1;
  }
  if(CheckDensity() == 2.7)
  {
    n1 = 1.37;
    n2 = 1;

  }
  if(VacuumToPlastic(r))
  {
    n1 =1;
    n2 = 1.58;
  }
  if(VacuumToAluminium(r))
  {
    n1 = 1;
    n2 = 1.37;


  }
  if(n1==n2){return false;};

  double theta = tools::Angle_Between_Vectors(Photons[i]->GetDirection(),n);
  if(CheckReflection(theta,n1,n2) )
  {
    vector<double> ndir = tools::Get_Reflected_Dir(Photons[i]->GetDirection(),n);
    geom->SetCurrentDirection(ndir.data());
    Photons[i]->ChangeDirection(ndir);

    return true;

  }else
  {
    vector<double> ndir = tools::Get_Refracted_Dir(Photons[i]->GetDirection(),n,theta,n1,n2);
    geom->SetCurrentDirection(ndir.data());
    Photons[i]->ChangeDirection(ndir);

    return false;
  }


}

bool Tracker::VacuumToPlastic(double r)
{
  if(CheckDensity()==0 && ((abs(r-Radius) < 1e-6) || (abs(abs(cpoint[2])-(0.5*Distance+Height))<1e-6) || (abs(abs(cpoint[2])-(0.5*Distance))<1e-6))){return true;};
  return false;
}

bool Tracker::VacuumToAluminium(double r)
{
  if(CheckDensity()==0 && ((abs(r-innerradius) < 1e-6) || (abs(r-outerradius) <1e-6) || (abs(abs(cpoint[2])-(Airgap+(0.5*Distance+Height)))<1e-6)|| (abs(abs(cpoint[2])-(-Airgap+(0.5*Distance)))<1e-6) || (abs(abs(cpoint[2])-(Thickness+Airgap+(0.5*Distance+Height)))<1e-6) || (abs(abs(cpoint[2])-(-Thickness-Airgap+(0.5*Distance)))<1e-6) ) ){return true;};
  return false;
}

//////////Muon Propagators//////////
//////////Muon Propagators//////////
//////////Muon Propagators//////////
//////////Muon Propagators//////////
//////////Muon Propagators//////////

void Tracker::Propagate_Muon()
{

  int scintillator_cross=0;

  while(!geom->IsOutside())
  {
    if(CheckDensity()==0)
    {
      Muon_Vacuum_Step();
    }
    if(CheckDensity()==1.023)
    {
      scintillator_cross++;
      Muon_Scintillator_Step();
    }
    if(CheckDensity()==2.7)
    {
      Muon_Aluminium_Step();
    }
    if(scintillator_cross == 2){DoubleCross=true;};

  }
  return;
}

void Tracker::Muon_Vacuum_Step()
{
  geom->FindNextBoundaryAndStep();
  Muon->IncreaseTime(geom->GetStep()/(Muon->GetVelocity()*2.998e10));
  Muon->ChangePosition(cpoint);
  geom->GetCurrentTrack()->AddPoint(cpoint[0],cpoint[1],cpoint[2],Muon->GetTime());
  return;
}

void Tracker::Muon_Scintillator_Step()
{
  geom->SetStep(stepvalue);
  while(CheckSameLocation())
  {
    geom->Step(kFALSE);
    Muon->IncreaseTime(stepvalue/(Muon->GetVelocity()*2.998e10));
    Muon->ChangePosition(cpoint);
    geom->GetCurrentTrack()->AddPoint(cpoint[0],cpoint[1],cpoint[2],Muon->GetTime());
    int n = generator->Generate_Photon_Number(10000*Update_Energy(stepvalue));
    N_photons+=n;
    for(int i=0;i<n;i++)
    {

      Photons.push_back(generator->Generate_Photon(Muon->GetPosition()));
      Photons[i]->IncreaseTime(Muon->GetTime());
    }
  }

  geom->FindNextBoundaryAndStep(stepvalue);
  Muon->IncreaseTime(geom->GetStep()/(Muon->GetVelocity()*2.998e10));
  Muon->ChangePosition(cpoint);
  geom->GetCurrentTrack()->AddPoint(cpoint[0],cpoint[1],cpoint[2],Muon->GetTime());
  int n = generator->Generate_Photon_Number(10000*Update_Energy(stepvalue));
  for(int i=0;i<n;i++)
  {
    //cout<<"Photon: "<<N_photons-n+i+1<<endl;
    Photons.push_back(generator->Generate_Photon(Muon->GetPosition()));
  }
  return;
}

void Tracker::Muon_Aluminium_Step()
{
  geom->FindNextBoundaryAndStep();
  Muon->IncreaseTime(1.37*geom->GetStep()/(Muon->GetVelocity()*2.998e10));
  Muon->ChangePosition(cpoint);
  geom->GetCurrentTrack()->AddPoint(cpoint[0],cpoint[1],cpoint[2],Muon->GetTime());
  return;
}

//////////Photon Propagators//////////
//////////Photon Propagators//////////
//////////Photon Propagators//////////
//////////Photon Propagators//////////
//////////Photon Propagators//////////

void Tracker::Propagate_Photons()
{
  for(int i=0;i<N_photons;i++)
  {

    InitializePhotonTrack(i);

    double absorption_step = generator->Generate_Photon_Step();
    double total_dist = 0;
    int j=0;
    while(!geom->IsOutside())
    {

      if(CheckDensity()==1.023)
      {
        geom->FindNextBoundary();
        if((total_dist+=geom->GetStep()) >=absorption_step)
        {
          Photon_Scintillator_Absorbtion(i,absorption_step-total_dist+geom->GetStep());
          break;
        }
        if(CheckDetector(cpoint))
        {
          N_detected++;
          break;
        }
        Photon_Scintillator_Step(i);

      }

      if(CheckDensity()==0)
      {
        geom->FindNextBoundary();
        Photon_Vacuum_Step(i);
      }

      if(CheckDensity()==2.7)
      {
        geom->FindNextBoundary();
        Photon_Aluminium_Step(i);
      }

    }

    if(geom->IsOutside()){N_lost++;};
  }

  return;
}

void Tracker::InitializePhotonTrack(int i)
{

  geom->SetCurrentTrack(geom->GetCurrentTrack()->AddDaughter(i,Photons[i]->GetPDG(),Photons[i]));
  geom->GetCurrentTrack()->AddPoint(Photons[i]->GetStartingPosition()[0],Photons[i]->GetStartingPosition()[1],Photons[i]->GetStartingPosition()[2],0);
  geom->InitTrack(geom->GetCurrentTrack()->GetFirstPoint(), Photons[i]->GetDirection().data());
  geom->SetCurrentPoint(Photons[i]->GetStartingPosition().data());
  return;
}

void Tracker::Photon_Scintillator_Step(int i)
{
  geom->Step(kTRUE,kFALSE);
  Photons[i]->IncreaseTime(1.58*geom->GetStep()/(2.998e10));
  Photons[i]->ChangePosition(cpoint);
  geom->GetCurrentTrack()->AddPoint(cpoint[0],cpoint[1],cpoint[2],Photons[i]->GetTime());

  if(!ReflectionHandler(i))
  {
    geom->FindNextBoundaryAndStep();
    Photons[i]->IncreaseTime(1.58*geom->GetStep()/(2.998e10));
    Photons[i]->ChangePosition(cpoint);
    geom->GetCurrentTrack()->AddPoint(cpoint[0],cpoint[1],cpoint[2],Photons[i]->GetTime());
  }
  return;
}

void Tracker::Photon_Scintillator_Absorbtion(int i, double step)
{
  geom->SetStep(step);
  geom->Step();
  Photons[i]->IncreaseTime(1.58*geom->GetStep()/2.998e10);
  Photons[i]->ChangePosition(cpoint);
  geom->GetCurrentTrack()->AddPoint(cpoint[0],cpoint[1],cpoint[2],Photons[i]->GetTime());
  N_absorbed++;
  return;
}

void Tracker::Photon_Vacuum_Step(int i)
{
  geom->Step(kTRUE,kFALSE);
  Photons[i]->IncreaseTime(geom->GetStep()/(2.998e10));
  Photons[i]->ChangePosition(cpoint);
  geom->GetCurrentTrack()->AddPoint(cpoint[0],cpoint[1],cpoint[2],Photons[i]->GetTime());

  if(!ReflectionHandler(i))
  {
    geom->FindNextBoundaryAndStep();
    Photons[i]->IncreaseTime(geom->GetStep()/(2.998e10));
    Photons[i]->ChangePosition(cpoint);
    geom->GetCurrentTrack()->AddPoint(cpoint[0],cpoint[1],cpoint[2],Photons[i]->GetTime());
  }
  return;
}

void Tracker::Photon_Aluminium_Step(int i)
{
  geom->Step(kTRUE,kFALSE);
  Photons[i]->IncreaseTime(1.37*geom->GetStep()/(2.998e10));
  Photons[i]->ChangePosition(cpoint);
  geom->GetCurrentTrack()->AddPoint(cpoint[0],cpoint[1],cpoint[2],Photons[i]->GetTime());

  if(!ReflectionHandler(i))
  {
    geom->FindNextBoundaryAndStep();
    Photons[i]->IncreaseTime(1.37*geom->GetStep()/(2.998e10));
    Photons[i]->ChangePosition(cpoint);
    geom->GetCurrentTrack()->AddPoint(cpoint[0],cpoint[1],cpoint[2],Photons[i]->GetTime());
  }
  return;
}

//////////Draw Mode/////////
//////////Draw Mode/////////
//////////Draw Mode/////////
//////////Draw Mode/////////
//////////Draw Mode/////////

void Tracker::Draw(){

    TCanvas *c1 = new TCanvas("c1","c1",1200,900);
    geom->GetTopVolume()->Draw();
    geom->DrawTracks("/*");
}

void Tracker::Propagate_Photons_DrawMode(int n)
{
  int m = N_photons/n;

  for(int j=0;j<n;j++)
  {
    int i=m*j;

    InitializePhotonTrack(i);
    double absorption_step = generator->Generate_Photon_Step();
    double total_dist = 0;

    while(!geom->IsOutside())
    {

      if(CheckDensity()==1.023)
      {
        geom->FindNextBoundary();
        if((total_dist+=geom->GetStep()) >=absorption_step)
        {
          Photon_Scintillator_Absorbtion(i,absorption_step-total_dist+geom->GetStep());
          break;
        }
        Photon_Scintillator_Step(i);
        if(CheckDetector(cpoint))
        {
          N_detected++;
          break;
        }
      }

      if(CheckDensity()==0)
      {
        geom->FindNextBoundary();
        Photon_Vacuum_Step(i);
      }

      if(CheckDensity()==2.7)
      {
        geom->FindNextBoundary();
        Photon_Aluminium_Step(i);
      }

    }

    if(geom->IsOutside()){N_lost++;};
  }
  cout<<endl<<"Total Photons Generated: "<<n<<endl;
  cout<<"Photons Absorbed: "<<N_absorbed<<endl;
  cout<<"Photons Detected: "<<N_detected<<endl;
  cout<<"Photons Lost: "<<N_lost<<endl;
  return;
}

/////////Debug Mode/////////
/////////Debug Mode/////////
/////////Debug Mode/////////
/////////Debug Mode/////////
/////////Debug Mode/////////

void Tracker::print_vector(const double* v)
{
  cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<endl;
}

void Tracker::Debug()
{
  N_photons = 1;
  for(int i=0;i<N_photons;i++)
  {
    cout<<"huh"<<endl;
    Photons.push_back(generator->Generate_Photon(generator->Generate_Position(Distance,Height-generator->Uniform(0,1)*Height,Radius)));
  }

  for(int i=0;i<N_photons;i++)
  {
    cout<<"Photon: "<<i<<endl;
    InitializePhotonTrack(i);

    double absorption_step = generator->Generate_Photon_Step();
    double total_dist = 0;
    cout<<"Absorbtion Step: "<<absorption_step<<endl;
    cout<<"Photon Energy: "<<Photons[i]->GetEnergy()<<endl<<endl;
    while(!geom->IsOutside())
    {

      if(CheckDensity()==1.023)
      {
        geom->FindNextBoundary();
        if((total_dist+=geom->GetStep()) >=absorption_step)
        {
          print_vector(geom->GetCurrentDirection());
          Photon_Scintillator_Absorbtion(i,absorption_step-total_dist+geom->GetStep());
          cout<<"Photon Absorbed"<<endl;
          break;
        }
        cout<<"Inside Scintillator"<<endl;
        cout<<"Distance traveled inside scintillators: "<<total_dist<<endl;
        cout<<"Position: "<<ends;
        print_vector(cpoint);
        cout<<"Direction: "<<ends;
        print_vector(geom->GetCurrentDirection());
        cout<<endl;

        Photon_Scintillator_Step(i);
        if(CheckDetector(cpoint))
        {
          print_vector(cpoint);
          cout<<"Photon Detected"<<endl;
          N_detected++;
          break;
        }

      }

      if(CheckDensity()==0)
      {
        cout<<"In Vacuum"<<endl;
        cout<<"Position: "<<ends;
        print_vector(cpoint);
        cout<<"Direction: "<<ends;
        print_vector(geom->GetCurrentDirection());
        cout<<endl;

        geom->FindNextBoundary();
        Photon_Vacuum_Step(i);
      }

      if(CheckDensity()==2.7)
      {
        cout<<"Inside Aluminium"<<endl;
        cout<<"Position: "<<ends;
        print_vector(cpoint);
        cout<<"Direction: "<<ends;
        print_vector(geom->GetCurrentDirection());
        cout<<endl;

        geom->FindNextBoundary();
        Photon_Aluminium_Step(i);
      }

    }
    cout<<"////////////////////"<<endl<<endl;
    if(geom->IsOutside()){N_lost++;};
  }

  cout<<endl<<"Total Photons Generated: "<<N_photons<<endl;
  cout<<"Photons Absorbed: "<<N_absorbed<<endl;
  cout<<"Photons Detected: "<<N_detected<<endl;
  cout<<"Photons Lost: "<<N_lost<<endl;
  Draw();
  return;

}
