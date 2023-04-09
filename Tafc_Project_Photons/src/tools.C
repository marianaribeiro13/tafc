#include "tools.h"
#include <numeric>

tools::~tools(){}


double tools::Norm(vector<double> v){
  double n = 0;
  for(int i=0; i < v.size(); i++){
    n+=v[i]*v[i];
  }

  return sqrt(n);
}

double tools::Angle_Between_Vectors(vector<double> v1,vector<double> v2){

  double dot=0;
  for(int i=0;i<3;i++)
  {
    dot+= v1[i]*v2[i];
  }
  double lenSq1 = v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
  double lenSq2 = v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2];

  return acos(dot/sqrt(lenSq1 * lenSq2));
}

double tools::SnellLaw(double thetai, double n1, double n2){

  return asin(n1*sin(thetai)/n2);
}

//Get reflected direction vector, from incident and normal vectors
vector<double> tools::Get_Reflected_Dir(vector<double> i, vector<double> n){

  std::vector<double> r(3);
  //double x = ;
  double x = 2*std::inner_product(i.begin(), i.end(), n.begin(), 0.);
  for(int j=0;j<3;j++){r[j] = (i[j]- x*n[j]);};
  //for(int j=0;j<3;j++){r[j]= -i[j];};
  return r;
}

//Get refracted direction vector, from incident and normal directions
vector<double> tools::Get_Refracted_Dir(vector<double> i, vector<double> n, double thetai, double n1, double n2)
{

  double r = n1 / n2;
  double c = std::inner_product(i.begin(), i.end(), n.begin(), 0.);
  double s = r*r*(1-c*c);
  vector<double> aux(3);
  for (int j=0; j<3; j++){
    aux[j] = r*i[j]+(r*c-sqrt(1-s))*n[j];
  }

  return tools::NormalizeVector(aux);
}

vector<double> tools::NormalizeVector(vector<double> v)
{
  double norm = Norm(v);
  for(int i=0;i<v.size();i++)
  {
    v[i]/=norm;
  }
  return v;
}


vector<string> tools::Read_File(string name){
  fstream fp;
  string aux;
  fp.open(name.c_str());
  if(fp.peek()!= EOF)
  {
    vector<string> v;
    while(!fp.eof())
    {
      getline(fp,aux);
      if( !(aux.find("//") == 0) && !(aux.find("#") == 0))
      {
        v.push_back(aux);
      }
    }
    return v;
  }
  fp.close();
  cout<<"File not Found"<<endl;
  exit(0);
}

TSpline3* tools::Interpolate_From_File(string name)
{
  vector<string> fs = tools::Read_File(name.c_str());
  vector<double> x,y;
  for(int i=0;i<fs.size()-1;i++)
  {
    double a,b;
    sscanf(fs[i].c_str(),"%lf %lf",&a ,&b);
    x.push_back(a);
    y.push_back(b);

  }
  auto I = new TSpline3("",x.data(),y.data(),x.size());
  return I;
}


void tools::print_vector(const double* v)
{
  cout<<v[0]<<" "<<v[1]<<" "<<v[2]<<endl;
  return;
}

double tools::RadialTheta(const double *cpoint)
{
  if(cpoint[0]>0)
  {
    return atan(cpoint[1]/cpoint[0]); 
  }else
  {
    if(cpoint[1]>=0) // x<0 and y>0
    {
      return atan(cpoint[1]/cpoint[0])+M_PI;
    }else // x<0 and y<0
    {
      return atan(cpoint[1]/cpoint[0])-M_PI;
    }
  }
  return 0;
}
