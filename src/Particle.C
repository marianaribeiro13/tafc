#include "Particle.h"
#include "TDatabasePDG.h"
#include <iostream>


////////////////////////////// Constructors ////////////////////////////////////

//If we know the momentum of the particle
Particle::Particle(int const PDG,  double p, std::vector<double> const& d) : pdg(PDG), momentum(p){

    //The TDatabasePDG particle database manager class creates a list of particles which by default is initialised from with the constants used by PYTHIA6
    //Get particle properties associated with given pdg
    TParticlePDG* part = TDatabasePDG::Instance()->GetParticle(pdg);

    //Get the particle mass in GeV
    double const m = part->Mass();

    //Convert mass to MeV
    mass = m*1000;

    //Compute energy in natural units
    energy = sqrt(momentum*momentum + mass*mass); //c=1

    direction = d;

}

//Copy Constructor
Particle::Particle(Particle* part) : pdg(part->pdg), mass(part->mass), energy(part->energy), momentum(part->momentum){

    direction = part->direction;
}

//////////////////////////////////////// Destructor //////////////////////////////////////////////

Particle::~Particle(){}

////////////////////////////////////// Getters /////////////////////////////////////////
std::vector<double> Particle::GetDirection(){

    std::vector<double> d {direction[0], direction[1], direction[2]};

    return d;
}
