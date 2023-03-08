#include <cmath>
#include <iostream>
#include <vector>
#include "Geometry.h"
#include "Muon.C"

using namespace std;

vector<double> generate_direction()
{
    vector<double> v(3);
    v[0] = (double)rand()/RAND_MAX -0.5;
    v[1] = (double)rand()/RAND_MAX -0.5;
    v[2] = abs((double)rand()/RAND_MAX -0.5);

    double norm = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    v[0] /= norm;
    v[1] /= norm;
    v[2] /= -norm;
    return v;
}

vector<double> generate_position(double R,double d,double h)
{
    vector<double> aux(3);
    aux[2] = d/2+h;
    double r = R * sqrt((double) rand()/ RAND_MAX);
    double theta = (double) rand()/RAND_MAX * 2 * 3.1415;
    aux[0] = r * cos(theta);
    aux[1] = r * sin(theta);
    return aux;

}





int main(int argc, char **argv)
{
  srand(time(0));
  double p[3]={0,0,0};
  double MM=0;
  double R = 25, h =1 , d=100;

  for(int i=0;i<1000;i++)
  {
      double e = 200*(double) rand()/RAND_MAX + 105.6583755;

      muon* M =  new muon(e,generate_position(R,d,h).data(),generate_direction().data());
      MM+=M->GetAngle()*(180/M_PI);

      for(int i=0;i<3;i++)
      {
          cout<<M->GetPosition()[i]<<" "<<ends;


      }

      cout<<endl<<M->GetAngle()*(180/M_PI)<<endl;

  }
  cout<<MM/1000<<endl;
}
