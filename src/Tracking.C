#include "Tracking.h"
#include "TCanvas.h"


////////////////////////// Constructor ////////////////////////////////////////

Tracking::Tracking(double radius, double height, double distance, double airgap, double althickness, double step, Generator* g) : Geometry(), stepvalue(step){
    
    generator = g;

    Build_MuonTelescope(radius, height, distance, airgap, althickness);
}


///////////////////////// Add Particle function //////////////////////////////

//Adds track associated to particle (give initial position)
int Tracking::AddParticle(int const id, vector<double> x, Particle* particle){

    std::cout << "Aquiii" << std::endl;

    //Create Track associated to particle
    int track_index = geom->AddTrack(id, particle->GetPDG(), particle); //track id, particle pdg, pointer to particle

    //Get created track
    TVirtualGeoTrack* track = geom->GetTrack(track_index);

    //Assign initial position to the track
    track->AddPoint(x[0], x[1], x[2], 0);

    return track_index;
}


/////////////////////////// Propagate Massive particle function //////////////////////////////

void Tracking::Propagate(int track_index){

    //Set the given track as the current track
    geom->SetCurrentTrack(track_index);

    //Get pointer to current track
    TVirtualGeoTrack* track = geom->GetCurrentTrack();

    //Get particle associated to track
    Particle* part = new Particle(dynamic_cast<Particle*>(track->GetParticle()));

    //Get particle direction
    std::vector<double> d = part->GetDirection();

    //Setting both initial point and direction and finding the state
    geom->InitTrack(track->GetFirstPoint(), d.data());

    //Propagation time
    double t = 0;

    while(!(geom->IsOutside())){

        //Get current material
        TGeoMaterial *cmat = CheckMaterial();

        //Check if the particle is in vacuum or not and propagate accordingly
        if(cmat->GetDensity() == 0){

            std::cout << "Outside scintilator" << endl;

            //Make step to the next boundary and cross it
            geom->FindNextBoundaryAndStep();

            //std::cout << part->GetEnergy() << std::endl;

            //Get the step taken
            double step_vac  = geom->GetStep();

            std::cout << "Step to cross boundary in vacuum: " << step_vac << std::endl;

            //Compute velocity
            double v = (double)(part->GetMomentum())/(part->GetEnergy());

            //Compute time
            t += double(step_vac)/v;

            //Get new position of the particle after making the step
            const double *cpoint = geom->GetCurrentPoint();

            //Add new point to the particle track
            track->AddPoint(*cpoint, *(cpoint+1), *(cpoint+2), t);

        } else {

            std::cout << "Inside scintilator" << endl;

            //Get current particle position in the geometry
            const double *cpoint = geom->GetCurrentPoint(); //cpoint always points to the current position

            //Number of arbitrary steps inside material
            int nsteps = 1;

            while(geom->IsSameLocation(*cpoint + stepvalue*d[0], *(cpoint+1) + stepvalue*d[1], *(cpoint+2) + stepvalue*d[2])){
                
                //Set arbitrary step
                geom->SetStep(stepvalue);

                //Execute step (flag = kFALSE means it is an arbitrary step, not limited by geometrical reasons)
                geom->Step(kFALSE);

                //Get energy before step
                double E = part->GetEnergy();

                //Compute particle velocity
                double v = (double)(part->GetMomentum())/E;

                //Compute propagation time
                t += double(stepvalue)/v;

                //Add new point to the particle track
                track->AddPoint(*cpoint, *(cpoint+1), *(cpoint+2),t);

                //Calculate energy lost to the material
                double dE = BetheBloch(v, stepvalue);

                //Update particle energy
                part->ChangeEnergy(E-dE);

                //std::cout << part->GetEnergy() << std::endl;

                //Update particle momentum
                double p = part->CalculateMomentum(E-dE);
                part->ChangeMomentum(p);

                //////////////////////////// Photons generation /////////////////////////////////

                //Get number of photons from Poisson distribution (approximatly 1 photon per 100 eV)
                int nphotons = generator->Generate_Photon_Number(dE*10000);

                //Save current position in vector
                vector<double> pos {*cpoint, *(cpoint+1), *(cpoint+2)};
                
                //if(nsteps == 100){

                std::cout << "Number of photons in step " << nsteps << ": " << nphotons << std::endl;

                for(int i = 0; i < nphotons; i++){

                    //int photon_index = AddParticle(i+1, pos, generator->Generate_Photon());

                    //generate photon according to Scintillator spectrum
                    Particle* photon = generator->Generate_Photon();

                    //Create secondary Track associated to photon
                    TVirtualGeoTrack* DaughterTrack = track->AddDaughter(i, photon->GetPDG(), photon); //track id, particle pdg, pointer to particle

                    //Assign initial position to the track
                    DaughterTrack->AddPoint(pos[0], pos[1], pos[2], t);
                
                    PropagatePhoton(DaughterTrack, t);
                }

                //}

                geom->SetCurrentTrack(track_index);

                geom->InitTrack(pos.data(), d.data());

                ++nsteps;
            }

            //Make step to the next boundary and cross it
            geom->FindNextBoundaryAndStep(stepvalue);

            //Get the step taken
            double step_mat  = geom->GetStep();

            //Get energy before step
            double E = part->GetEnergy();

            //Compute velocity
            double v = (double)(part->GetMomentum())/E;

            //Compute time
            t += double(step_mat)/v;

            //Add new point to the particle track
            track->AddPoint(*cpoint, *(cpoint+1), *(cpoint+2), t);

            //Calculate energy lost to the material
            double dE = BetheBloch(v, step_mat);

            std::cout << "Step to cross boundary inside scintillator: " << step_mat << std::endl;

            //Update particle energy
            part->ChangeEnergy(E-dE);

            //std::cout << part->GetEnergy() << std::endl;

            //Update particle momentum
            double p = part->CalculateMomentum(E-dE);
            part->ChangeMomentum(p);

        }
    }
}


//////////////////////////////////// Propagate Photon function ///////////////////////////////

void Tracking::PropagatePhoton(TVirtualGeoTrack* track, double t){

    //Set the given track as the current track
    geom->SetCurrentTrack(track);

    //Get particle associated to track
    Particle* part = new Particle(dynamic_cast<Particle*>(track->GetParticle()));

    //Get particle direction
    std::vector<double> d = part->GetDirection();

    //Setting both initial point and direction and finding the state
    geom->InitTrack(track->GetFirstPoint(), d.data());

    //Get new position of the particle after making the step
    const double *cpoint = geom->GetCurrentPoint();

    while(!(geom->IsOutside())){

        //Get current material
        TGeoMaterial *cmat = CheckMaterial();

        //Check if the particle is in vacuum or not and propagate accordingly
        if(cmat->GetDensity() == 0){

            std::cout << "Outside scintilator" << endl;

            //Make step to the next boundary
            geom->FindNextBoundary();

            //Get the step taken
            double step_vac  = geom->GetStep();

            if(Is_Reflected(d)){
                //Dont cross boundary
                geom->Step(kTRUE, kFALSE);
                
                //new direction

            } else {
                //Cross boundary
                geom->Step(kTRUE, kTRUE);

                //newdirection
            }

            //Compute time
            t += double(step_vac); //v=c=1

            //Add new point to the particle track
            track->AddPoint(*cpoint, *(cpoint+1), *(cpoint+2), t);

        } else { 

            //Generate random step according to the probability of absorption of the photon
            double absorption_step = generator->Generate_Photon_Step();

            double total_dist = 0.;

            while(total_dist < absorption_step){
                //if photon step 
                geom->FindNextBoundary();

                double snext = geom->GetStep();

                if(absorption_step < snext){

                    //Set arbitrary step
                    geom->SetStep(absorption_step);

                    //Execute step (flag = kFALSE means it is an arbitrary step, not limited by geometrical reasons)
                    geom->Step(kFALSE);

                    //Compute time
                    t += double(snext)/v;

                    //Add new point to the particle track
                    track->AddPoint(*cpoint, *(cpoint+1), *(cpoint+2), t);

                    return;

                } else {
                    if(Is_Reflected(d)){
                        //Dont cross boundary
                        geom->Step(kTRUE, kFALSE);
                        
                        //new direction
                        total_dist += geom->GetStep();

                        //Compute time
                        t += double(snext)/v;

                        //Add new point to the particle track
                        track->AddPoint(*cpoint, *(cpoint+1), *(cpoint+2), t);

                    } else {
                        //Cross boundary
                        geom->Step(kTRUE, kTRUE);
                        
                        //Get refraction angle
                        double thetat = SnellLaw(thetai, 1, 1.58);

                        //newdirection
                        

                        //Compute time
                        t += double(snext)/v;

                        //Add new point to the particle track
                        track->AddPoint(*cpoint, *(cpoint+1), *(cpoint+2), t);

                        break;
                    }
                }
            }
        }
    }
}

 //double thetat = SnellLaw(thetai, 1, 1.58);

        //new direction;

        /*//Make step to the next boundary and cross it
        geom->FindNextBoundaryAndStep();

        //Get the step taken
        double snext  = geom->GetStep();

        //std::cout << "Photon Step to cross boundary in vacuum: " << snext << std::endl;

        
   // }
//}


//////////////////////////////////// Check current material /////////////////////////////

TGeoMaterial* Tracking::CheckMaterial(){

    //Get current node
    TGeoNode *cnode = geom->GetCurrentNode();

    //Get current volume
    TGeoVolume *cvol = cnode->GetVolume();

    //Get current material
    TGeoMaterial *cmat = cvol->GetMedium()->GetMaterial();

    return cmat;
}


////////////////////////////////////// Bethe Bloch function to calculate energy loss in material //////////////////////////////////

double Tracking::BetheBloch(double v, double step){

    //NATURAL
    //Double_t alpha=0.0072974;
    //Double_t I1=64.70713358e-6; //MeV
    //Double_t n=2.5587248e18; //MeV
    Int_t Z=-1;
    //Double_t mass_muon=105.6583755; //MeV
    //Double_t mass_electron=0.51099895; //MeV
    //Double_t dE=0;
    //Double_t dx=1.e-4;

    //S.I.
    Double_t c = 299792458;
    Double_t mp = 1.672621637e-27;
    Double_t me = 9.1093821499999992e-31;
    Double_t qe = 1.602176487e-19;
    Double_t na = 6.02214179e23;
    Double_t eps0 = 8.854187817e-12;
    Double_t n_density=3.33e29; // per m^3
    Double_t I = 1.03660828e-17;

    //double dEdx = (4*TMath::Pi()*n*Z*Z*alpha*alpha*(log((2*mass_electron*v*v)/(I*(1.-v*v)))-v*v))/(mass_electron*v*v);

    double dEdx_SI = (qe*qe*qe*qe*n_density*Z*Z*(log((2*me*c*c*v*v)/(I*(1-v*v)))-v*v))/(4*M_PI*eps0*eps0*me*c*c*v*v);
    dEdx_SI = dEdx_SI/(1.602e-13);

    double dE_SI = dEdx_SI * step/100;

    return dE_SI;
}


double Tracking::FresnelLaw(double thetai, double n1, double n2){
    
    // Reflection probability for s-polarized light
    double Rs = abs((n1*cos(thetai)-n2*sqrt(1-(n1*sin(thetai)/n2)*(n1*sin(thetai)/n2)))/
                    (n1*cos(thetai)+n2*sqrt(1-(n1*sin(thetai)/n2)*(n1*sin(thetai)/n2))));
    
    // Reflection probability for p-polarized light
    double Rp = abs((n1*sqrt(1-(n1*sin(thetai)/n2)*(n1*sin(thetai)/n2))-n2*cos(thetai)))/
                    (n1*sqrt(1-(n1*sin(thetai)/n2)*(n1*sin(thetai)/n2))+n2*cos(thetai));

    return 0.5*(Rs+Rp);
}

//Check if light is reflected or transmitted and get new light direction
bool Tracking::Check_Reflection(std::vector<double>& di){

    double* normal = geom->FindNormal(Bool_t forward=kTRUE);

    vector<double> n {*normal, *(normal+1), *(normal+2)};

    double thetai = tools::Angle_Between_Vectors(di, n);

    double Reff = FresnelLaw(thetai, 1, 1.58);

    if(generator->Uniform(0,1) < Reff) {

        return true;

    } else {

        return false;
    }
}


////////////////////////////////////////////////// Draw geometry and tracks function //////////////////////////////////////

void Tracking::Draw(){

    TCanvas *c1 = new TCanvas("c1","c1");
    //c1->cd();
    //top->Draw("ogle");
    geom->GetTopVolume()->Draw();
    geom->DrawTracks("/*");
    //geom->AnimateTracks(0, 1000, 200, "/* /G /S");
    c1->SaveAs("Simulation.pdf");
}
