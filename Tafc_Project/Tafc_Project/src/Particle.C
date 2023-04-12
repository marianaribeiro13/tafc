#include "Particle.h"
#include "TDatabasePDG.h"
#include <iostream>


////////////////////////////// Constructors ////////////////////////////////////

//If we know the momentum of the particle
Particle::Particle(int const PDG,  double p, std::vector<double> const& d,vector<double> const& x) : pdg(PDG), momentum(p){

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
    StartingPosition = x;
    Position = x;
    time = 0;

}

//If we know the momentum of the particle
Particle::Particle(int const PDG,  double m, double p, std::vector<double> const& d,vector<double> const& x) : pdg(PDG), momentum(p){
    
    //Convert mass to MeV
    mass = m*1000;

    //Compute energy in natural units
    energy = sqrt(momentum*momentum + mass*mass); //c=1

    direction = d;
    StartingPosition = x;
    Position = x;
    time = 0;

}

//Copy Constructor
Particle::Particle(Particle* part) : pdg(part->pdg), mass(part->mass), energy(part->energy), momentum(part->momentum){

    direction = part->direction;
}

void Particle::ChangePosition(const double *cpoint)
{
    for(int i=0;i<3;i++)
    {
        Position[i] = cpoint[i];
    }
    return;
}

//////////////////////////////////////// Destructor //////////////////////////////////////////////

Particle::~Particle(){}
