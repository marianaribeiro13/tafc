#include "tools.h"

tools::~tools(){}

vector<double> tools::Generate_Vector()
{
  vector<double> v(3);
  v[0] = (double)rand()/RAND_MAX -0.5;
  v[1] = (double)rand()/RAND_MAX -0.5;
  v[2] = (double)rand()/RAND_MAX -0.5;

  double norm = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  v[0] /= norm;
  v[1] /= norm;
  v[2] /= norm;
  return v;
}

double tools::Random_Distribution(double xmin,double xmax,TF1 *F)
{
  double x = xmin+(xmax-xmin)*(double) rand()/RAND_MAX;
  double y = F->GetMaximum(xmin,xmax)*(double) rand()/RAND_MAX;

  while(y>F->Eval(x))
  {
    x = xmin+(xmax-xmin)*(double) rand()/RAND_MAX;
    y = F->GetMaximum(xmin,xmax)*(double) rand()/RAND_MAX;

  }

  return x;
}

muon* tools::Generate_Muon(vector<double> Position)
{
  vector<double> Direction(3);
  double angle = 0.5*M_PI * (double) rand()/RAND_MAX;
  double aux0 = (double) rand()/RAND_MAX,aux1 = (double) rand()/RAND_MAX;
  Direction[2] = -cos(angle);
  Direction[0] = aux0*sqrt(1-Direction[2]*Direction[2])/(sqrt(aux0*aux0+aux1*aux1));
  Direction[1] = aux1*sqrt(1-Direction[2]*Direction[2])/(sqrt(aux0*aux0+aux1*aux1));
  double energy = 1000;
  //auto f = [](double *x double *par){return pow(cos(x[1]),3)*0.00253*pow(x[0]*cos(x[1]),-(0.2455+1.288*log10(x[0]*x[1])));};

  cout<<Direction[0]<<" "<<Direction[1]<<" "<<Direction[2]<<" "<<tools::Norm(Direction)<<endl;

  muon *M = new muon(energy,Position.data(),Direction.data());

  return M;
}

double tools::Norm(vector<double> v)
{
  double n = 0;
  for(int i=0;i<v.size();i++)
  {
    n+=v[i]*v[i];
  }
  return sqrt(n);
}

vector<double> tools::Generate_Position(double R,double d,double h)
{
    vector<double> aux(3);
    aux[2] = d/2+h;
    double r = R * sqrt((double) rand()/ RAND_MAX);
    double theta = (double) rand()/RAND_MAX * 2 * 3.1415;
    aux[0] = r * cos(theta);
    aux[1] = r * sin(theta);
    return aux;

}

vector<double> tools::Random_Distribution_2D(TF1* F,double xmin,double xmax,double ymin,double ymax,double max)
{
    vector<double> x(2);
    x[0]= xmin+(xmax-xmin)*(double) rand()/RAND_MAX;
    x[1]= ymin+(ymax-ymin)*(double) rand()/RAND_MAX;
    double y = max*(double) rand()/RAND_MAX;
    while(y>F->EvalPar(x.data()))
    {
        x[0]= xmin+(xmax-xmin)*(double) rand()/RAND_MAX;
        x[1]= ymin+(ymax-ymin)*(double) rand()/RAND_MAX;
        y = max*(double) rand()/RAND_MAX;
    }
    return x;
}
