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

double tools::Angle_Between_Vectors(vector<double>& v1, vector<double>& v2){
  
  double dot = std::inner_product(v1.begin(), v1.end(), v2.begin(), 0);   
  double lenSq1 = v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
  double lenSq2 = v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2];
  
  return acos(dot/sqrt(lenSq1 * lenSq2));
}

double tools::SnellLaw(double thetai, double n1, double n2){

  return asin(n1*sin(thetai)/n2);
}

//Get reflected direction vector, from incident and normal vectors
vector<double> tools::Get_Reflected_Dir(vector<double>& i, vector<double>& n){

  std::vector<double> r(3);
  double x = 2*std::inner_product(i.begin(), i.end(), n.begin(), 0);
  for(int j=0;j<3;j++){r[j]=i[j]-x*n[j];};
  return r;
}

//Get refracted direction vector, from incident and normal directions
vector<double> tools::Get_Refracted_Dir(vector<double>& i, vector<double>& n, double thetai, double n1, double n2){
  
  const double r = n1 / n2;
  const double sinT2 = r * r * (1.0 - thetai * thetai);
  if(sinT2 > 1.0) exit(0); // TIR
  const double cosT = sqrt(1.0 - sinT2);
  vector<double> aux(3);
  for (i=1; i<3; i++){
    aux[i]=r*i[i]+(-r*thetai -cosT)*n[i];
  }
  return aux;
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

TSpline3* tools::Interpolate_Photon_Spectrum(string name){
  vector<string> fs = tools::Read_File(name.c_str());
  vector<double> x,y;
  for(int i=0;i<fs.size();i++)
  {
    double a,b;
    sscanf(fs[i].c_str(),"%lf %lf",&a ,&b);
    x.push_back(a);
    y.push_back(b);
  }
  auto I = new TSpline3("f",x.data(),y.data(),x.size());
  return I;
}
