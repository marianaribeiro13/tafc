#include "Tracking.h"
#include "TApplication.h"

int main(){

    Generator* gen = new Generator();

    double distance = 25.;
    double step = 0.001;
    
    //Define geometry and tracking step
    Tracking Simulation(distance, step, gen);

    //std::cout << "Aquiii" << std::endl;

    //Generate random incident muon
    Particle* muon = Simulation.GenerateCosmicMuon();

    //Particle* muon2 = Simulation.GenerateCosmicMuon();

    // Generate random initial position of the muon (uniform distribution along the incident plane)
    std::vector<double> x = gen->Generate_Position(distance);

    //std::vector<double> x1 = gen->Generate_Position(distance);

    //Add muon to the generated position and create associated track
    int muonTrack_index = Simulation.AddParticle(0, x, muon);

    //int muonTrack_index2 = Simulation.AddParticle(1, x1, muon2);

    //Propagate muon

    Simulation.Propagate(muonTrack_index);

    //Simulation.Propagate(muonTrack_index2);

    TApplication app("app", nullptr, nullptr);

    //Draw telescope and tracks
    Simulation.Draw();

    app.Run();

    return 0;
}