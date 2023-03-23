#include "Generators.h"

Generator::Generator()
{
   Random = new TRandom(time(0));
   auto f = [](double *x,double *par)
   {
       return pow(cos(x[1]),3)*0.00253*pow(x[0]*cos(x[1]),-(0.2455+1.288*log10(x[0]*cos(x[1]))-0.255*pow(log10(x[0]*cos(x[1])),2)+0.0209*pow(log10(x[0]*cos(x[1])),3) ));
   };
   Momentum_Distribution = new TF1("f",f);
   Photon_Spectrum = tools::Interpolate_Photon_Spectrum("Photon_Spectrum.txt");
}

Generator::~Generator(){};

vector<double> Generator::Generate_Vector()
{
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

double Generator::Random_Distribution(double xmin,double xmax,TF1 *F)
{
  double x = Random->Uniform(xmin,xmax);
  double y = F->GetMaximum(xmin,xmax)*Random->Uniform(1);

  while(y>F->Eval(x))
  {
     x = Random->Uniform(xmin,xmax);
     y = F->GetMaximum(xmin,xmax)*Random->Uniform(1);

  }

  return x;
}

vector<double> Generator::Generate_Position(double d)
{
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
    while(y>F->EvalPar(x.data()))
    {
      x[0]= Random->Uniform(xmin,xmax);
      x[1]= Random->Uniform(ymin,ymax);
      y = Random->Uniform(max);
    }
    return x;
}

double Generator::Generate_Photon_Energy()
{
  double x = Random->Uniform(380,500);
  double y = Random->Uniform(1);

  while(y>Photon_Spectrum->Eval(x))
  {
    x = Random->Uniform(380,500);
    y = Random->Uniform(1);
  }
  return (6.626e-34*2.998e8)/(x*1e-9)/1.602e-19; //Convert to eV
}

int Generator::Generate_Photon_Number(double expected)
{
  return Random->Poisson(expected);
}

// Muon* Generator::Generate_Muon(double distance)
// {
//   vector<double> aux = Random_Distribution_2D(Momentum_Distribution,1,2000,0,M_PI/2.,Momentum_Distribution->GetMaximum());
//   vector<double> direction(3);
//   direction[2] = -cos(aux[1]);
//   double aux0 = Random->Uniform(-1,1), aux1 = Random->Uniform(-1,1);
//   direction[0] = aux0*sqrt((1-direction[2]*direction[2])/(aux0*aux0+aux1*aux1));
//   direction[1] = aux1*sqrt((1-direction[2]*direction[2])/(aux0*aux0+aux1*aux1));
//   Muon *M = new Muon(Generate_Position(distance),direction,1000*aux[0]);
//   return M;
// }
//
// Photon* Generator::Generate_Photon(vector<double> x)
// {
//
//   return new Photon(x,Generate_Vector(),Generate_Photon_Energy());
// }
