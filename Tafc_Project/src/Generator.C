#include "Generator.h"
#include "TMath.h"
#include <iostream>

Generator::Generator(){
  Random = new TRandom(time(0));

   //2D Spectrum of the muon (x[0] is the momentum and x[1] is the incident angle)
  auto f = [](double *x,double *par){
       return pow(cos(x[1]),3)*0.00253*pow(x[0]*cos(x[1]),-(0.2455+1.288*log10(x[0]*cos(x[1]))-0.255*pow(log10(x[0]*cos(x[1])),2)+0.0209*pow(log10(x[0]*cos(x[1])),3) ));
  };

  //Create TF1 function with the 2D spectrum of the muon
  Momentum_Distribution = new TF1("f",f);
  Photon_Spectrum = tools::Interpolate_From_File("Photon_Spectrum.txt");
  Detector_Efficiency = tools::Interpolate_From_File("Detector_Efficiency_Spectrum.txt");
}

Generator::Generator(int seed)
{
 
  Random = new TRandom(seed);

   //2D Spectrum of the muon (x[0] is the momentum and x[1] is the incident angle)
  auto f = [](double *x,double *par){
       return pow(cos(x[1]),3)*0.00253*pow(x[0]*cos(x[1]),-(0.2455+1.288*log10(x[0]*cos(x[1]))-0.255*pow(log10(x[0]*cos(x[1])),2)+0.0209*pow(log10(x[0]*cos(x[1])),3) ));
  };

  //Create TF1 function with the 2D spectrum of the muon
  Momentum_Distribution = new TF1("f",f);

  Photon_Spectrum = tools::Interpolate_From_File("Photon_Spectrum.txt");
  Detector_Efficiency = tools::Interpolate_From_File("Detector_Efficiency_Spectrum.txt");
}


Generator::~Generator(){};

vector<double> Generator::Generate_Vector(){

  vector<double> v(3);
  v[0] = Random->Uniform(-1,1);
  v[1] = Random->Uniform(-1,1);
  v[2] = Random->Uniform(-1,1);

  double norm = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  v[0] /= norm;
  v[1] /= norm;
  v[2] /= norm;
  return v;

}

vector<double> Generator::Generate_Direction_From_Theta(double theta) {

  vector<double> d(3);

  d[2] = -cos(theta); // z direction according to theta

  //generate 2 random numbers between -1 and 1
  double aux0 = Uniform(-1,1), aux1= Uniform(-1,1);

  //generate x and y direction normalized
  d[0] = aux0*(double)sqrt((1-d[2]*d[2])/(aux0*aux0+aux1*aux1)); // x direction
  d[1] = aux1*(double)sqrt((1-d[2]*d[2])/(aux0*aux0+aux1*aux1)); // y direction

  return d;
}

double Generator::Random_Distribution(double xmin,double xmax,TF1 *F){

  double x = Random->Uniform(xmin,xmax);
  double y = F->GetMaximum(xmin,xmax)*Random->Uniform(1);

  while(y>F->Eval(x))
  {
     x = Random->Uniform(xmin,xmax);
     y = F->GetMaximum(xmin,xmax)*Random->Uniform(1);
  }

  return x;

}

vector<double> Generator::Generate_Position(double d,double h,double R){


  vector<double> aux(3);
  aux[2] = d/2+h-1e-12;
  double r = Random->Uniform(0,R);
  double theta = Random->Uniform(0,2*M_PI);
  aux[0] = r * cos(theta);
  aux[1] = r * sin(theta);
  return aux;

}

vector<double> Generator::Random_Distribution_2D(TF1* F,double xmin,double xmax,double ymin,double ymax,double max)
{
    vector<double> x(2);
    x[0]= Random->Uniform(xmin,xmax);
    x[1]= Random->Uniform(ymin,ymax);
    double y = Random->Uniform(max);

    while(y > F->EvalPar(x.data())) {
      x[0]= Random->Uniform(xmin,xmax);
      x[1]= Random->Uniform(ymin,ymax);
      y = Random->Uniform(max);
    }

    return x;
}

double Generator::Generate_Photon_Energy(){

  double x = Random->Uniform(380,500);
  double y = Random->Uniform(1);

while(y>Photon_Spectrum->Eval(x)){

    x = Random->Uniform(380,500);
    y = Random->Uniform(1);

  }

  return (6.626e-34*2.998e8)/(x*1e-9)/1.602e-13; //Convert to MeV
}

int Generator::Generate_Photon_Number(double expected){

  return Random->Poisson(expected);
}

Particle* Generator::Generate_CosmicMuon(vector<double> x)
{

  //Generate random momentum and incident angle according to 2D pdf using acceptance-rejection
  vector<double> aux = Random_Distribution_2D(GetMomentumDistribution(),1,2000,0, TMath::Pi()/2, GetMomentumDistribution()->GetMaximum());

  //Create muon (pdg = 13)
  Particle* muon = new Particle(13, 0.105658, aux[0]*1000, Generate_Direction_From_Theta(aux[1]),x);

  return muon;
}

Particle* Generator::Generate_Photon(vector<double> x){

  return new Particle(22 , 0, Generate_Photon_Energy(), Generate_Vector(),x);
}



double Generator::Generate_Photon_Step(){

  return -380*log(Random->Uniform(1));
}
