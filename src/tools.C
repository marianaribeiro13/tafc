#include "tools.h"

tools::~tools(){}


double tools::Norm(vector<double> v)
{
  double n = 0;
  for(int i=0;i<v.size();i++)
  {
    n+=v[i]*v[i];
  }
  return sqrt(n);
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
