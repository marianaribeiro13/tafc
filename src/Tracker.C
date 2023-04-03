#include "Tracker.h"

Tracker::Tracker(double radius, double height, double distance, double airgap, double althickness, double step, Generator* g, int n_SIPMS,double SIPM_size)
: Geometry(), stepvalue(step), N_photons(0),N_absorbed(0),N_detected(0),N_lost(0),DoubleCross(false)
{
  generator = g;
  Build_MuonTelescope(radius, height, distance, airgap, althickness,n_SIPMS,SIPM_size);
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

Tracker::~Tracker()
{

  delete Muon;
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
  if((abs(h-(0.5*Distance+Height))<1e-6) || (abs(h-0.5*Distance)<1e-6) || (abs(h-(Airgap+(0.5*Distance)+Height))<1e-6)  || (abs(h-(-Airgap+(0.5*Distance)))<1e-6) || (abs(h-(Thickness+Airgap+(0.5*Distance)+Height))<1e-6) || (abs(h-(-Thickness-Airgap+(0.5*Distance)))<1e-6)  )
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

bool Tracker::VacuumToPlastic(double r)
{
  if(CheckDensity()==0 && ((abs(r-Radius) < 1e-6) || (abs(abs(cpoint[2])-(0.5*Distance+Height))<1e-6) || (abs(abs(cpoint[2])-(0.5*Distance))<1e-6))){return true;};
  return false;
}

bool Tracker::VacuumToAluminium(double r)
{
  if(CheckDensity()==0 && ((abs(r-innerradius) < 1e-6) || (abs(r-outerradius) <1e-6) || (abs(abs(cpoint[2])-(Airgap+0.5*Distance+Height))<1e-6)|| (abs(abs(cpoint[2])-(-Airgap+(0.5*Distance)))<1e-6) || (abs(abs(cpoint[2])-(Thickness+Airgap+0.5*Distance+Height))<1e-6) || (abs(abs(cpoint[2])-(-Thickness-Airgap+(0.5*Distance)))<1e-6) ) ){return true;};
  return false;
}

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

    Photons.push_back(generator->Generate_Photon(Muon->GetPosition()));
  }
  return;
}

void Tracker::Muon_Aluminium_Step()
{
  geom->FindNextBoundaryAndStep();
  Muon->IncreaseTime(geom->GetStep()/(Muon->GetVelocity()*2.998e10));
  Muon->ChangePosition(cpoint);
  geom->GetCurrentTrack()->AddPoint(cpoint[0],cpoint[1],cpoint[2],Muon->GetTime());
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

    double absorption_step = generator->Generate_Photon_Step();
    double total_dist = 0;
    int l =0;
    while(true)
    {

      geom->FindNextBoundary();
      geom->Step(kTRUE,kFALSE);

      if(CheckDensity()==1.023)
      {

        if((total_dist+=geom->GetStep()) >=absorption_step)
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

      Photons[i]->IncreaseTime(GetRefractiveIndex()*geom->GetStep()/(2.998e10));
      Photons[i]->ChangePosition(cpoint);
      geom->GetCurrentTrack()->AddPoint(cpoint[0],cpoint[1],cpoint[2],Photons[i]->GetTime());
    }
    delete Photons[i];
    //if(geom->IsOutside()){N_lost++;};

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

void Tracker::Photon_Scintillator_Reflection_Check(int i)
{
  vector<double> n = GetNormal();
  double theta = tools::Angle_Between_Vectors(Photons[i]->GetDirection(),n);
  if(CheckReflection(theta,1.58,1) )
  {
    vector<double> ndir = tools::Get_Reflected_Dir(Photons[i]->GetDirection(),n);
    geom->SetCurrentDirection(ndir.data());
    Photons[i]->ChangeDirection(ndir);

  }else
  {
    vector<double> ndir = tools::Get_Refracted_Dir(Photons[i]->GetDirection(),n,theta,1.58,1);
    geom->FindNextBoundaryAndStep();
    geom->SetCurrentDirection(ndir.data());
    Photons[i]->ChangeDirection(ndir);

  }
  return;
}

bool Tracker::Photon_Vacuum_Reflection_Check(int i)
{

  vector<double> n = GetNormal();
  double theta = tools::Angle_Between_Vectors(Photons[i]->GetDirection(),n);
  double r = sqrt(cpoint[0]*cpoint[0]+cpoint[1]*cpoint[1]);
  if(VacuumToPlastic(r))
  {
    if(CheckReflection(theta,1,1.58) )
    {
      vector<double> ndir = tools::Get_Reflected_Dir(Photons[i]->GetDirection(),n);
      geom->SetCurrentDirection(ndir.data());
      Photons[i]->ChangeDirection(ndir);
    }else
    {
      vector<double> ndir = tools::Get_Refracted_Dir(Photons[i]->GetDirection(),n,theta,1,1.58);
      geom->SetCurrentDirection(ndir.data());
      Photons[i]->ChangeDirection(ndir);
      geom->FindNextBoundaryAndStep();
    }
  }else
  {
    if(VacuumToAluminium(r))
    {
      if(generator->Uniform(0,1) < .92)
      {
        vector<double> ndir = tools::Get_Reflected_Dir(Photons[i]->GetDirection(),n);
        geom->SetCurrentDirection(ndir.data());
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