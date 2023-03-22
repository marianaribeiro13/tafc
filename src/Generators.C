#include "Generators.h"

Generator::Generator()
{
   Random = new TRandom(time(0));
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

vector<double> Generator::Generate_Position(double R,double d,double h)
{
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

double Generator::Generate_Photon_Energy(TSpline3* F)
{
  double x = Random->Uniform(380,500);
  double y = Random->Uniform(1);

  while(y>F->Eval(x))
  {
    x = Random->Uniform(380,500);
    y = Random->Uniform(1);
  }
  return x;
}

int Generator::Generate_Photon_Number(double expected)
{
  return Random->Poisson(expected);
}
