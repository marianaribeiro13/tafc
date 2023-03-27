#include "Generators.h"
#include "TMath.h"

Generator::Generator(){
  Random = new TRandom(time(0));

   //2D Spectrum of the muon (x[0] is the momentum and x[1] is the zenith angle)
  auto f = [](double *x,double *par){
       return pow(cos(x[1]),3)*0.00253*pow(x[0]*cos(x[1]),-(0.2455+1.288*log10(x[0]*cos(x[1]))-0.255*pow(log10(x[0]*cos(x[1])),2)+0.0209*pow(log10(x[0]*cos(x[1])),3) ));
  };

  //Create TF1 function with the 2D spectrum of the muon
  Momentum_Distribution = new TF1("f",f);

  //Interpolate the data obtained for the photon spectrum
  Photon_Spectrum = tools::Interpolate_Photon_Spectrum("Photon_Spectrum.txt");

  auto f2 = [](double *x,double *par)
  {
    return exp(-x[0]/380.);
  };
  Absorbtion_Probability = new TF1("f2",f2);
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

vector<double> Generator::Generate_Position(double d){

  double h=1,R=5.;
  vector<double> aux(3);
  aux[2] = d/2+h;
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

Particle* Generator::Generate_CosmicMuon(){

  //Generate random momentum and incident angle according to 2D pdf using acceptance-rejection
  vector<double> aux = Random_Distribution_2D(GetMomentumDistribution(),1,2000,0, TMath::Pi()/2, GetMomentumDistribution()->GetMaximum());

  //Create muon (pdg = 13)
  Particle* muon = new Particle(13, aux[0]*1000, Generate_Direction_From_Theta(aux[1]));

  return muon;
}

Particle* Generator::Generate_Photon(){

  return new Particle(22 ,Generate_Photon_Energy(), Generate_Vector());
}

  double Generator::Generate_Photon_Step()
{
  return 380*log(1/(1-Random->Uniform(1)));

}
