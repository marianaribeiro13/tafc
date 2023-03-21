#include "Muon.h"

//Considering only muons in the x-z plane
Muon::Muon(std::vector<double> x, double theta, double E) : pdg(13), energy(E){
  
  momentum = sqrt(energy*energy - mass*mass); //c=1
  velocity = double(momentum)/energy; //c=1

  direction.reserve(3);
  direction[0] = sin(theta);
  direction[1] = 0;
  direction[2] = -cos(theta);

  position.reserve(4);
  for(int i=0; i<4; i++) {
    position[i]= x[i];
  }

}

Muon::~Muon(){}

std::vector<double> Muon::GetPosition(){

  std::vector<double> x {position[0], position[1], position[2], position[3]};

  return x;
}

std::vector<double> Muon::GetDirection(){

  std::vector<double> d {direction[0], direction[1], direction[2]};

  return d;
}
