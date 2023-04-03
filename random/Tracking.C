#include "Tracking.h"
#include "TCanvas.h"


////////////////////////// Constructor ////////////////////////////////////////

Tracking::Tracking(double radius, double height, double distance, double airgap, double althickness, double step, Generator* g) : Geometry(), stepvalue(step){

    generator = g;

    Build_MuonTelescope(radius, height, distance, airgap, althickness);
}


///////////////////////// Add Particle function //////////////////////////////

//Adds track associated to particle (give initial position)
//id is a track identifier that can be useful
int Tracking::AddParticle(int const id, vector<double> x, Particle* particle){

    //Create Track associated to particle object
    int track_index = geom->AddTrack(id, particle->GetPDG(), particle); //track id, particle pdg, pointer to particle

    //Get created track
    TVirtualGeoTrack* track = geom->GetTrack(track_index);

    //Assign initial position and time to the track
    track->AddPoint(x[0], x[1], x[2], 0);

    return track_index;
}


/////////////////////////// Propagate Massive particle function //////////////////////////////

//This propagation function considers particles which (on average) do not change direction when interacting with the medium
//The energy loss of the particle through the material is calculated through the Bethe-Bloch equation
void Tracking::Propagate(int track_index){

    //Set the given track as the current track
    geom->SetCurrentTrack(track_index);

    //Get pointer to current track
    TVirtualGeoTrack* track = geom->GetCurrentTrack();

    //Get particle associated to track
    Particle* part = new Particle(dynamic_cast<Particle*>(track->GetParticle()));

    //Get particle current direction
    std::vector<double> d = part->GetDirection();

    //Setting both initial point and direction and finding the state
    geom->InitTrack(track->GetFirstPoint(), d.data());

    //Propagation time
    double t = 0;

    //Velocity variable
    double v = 0;

    //Get pointer to current position (Whenever we update the position cpoint is updated,
    //since it is equal to the geom pointer to the current position)
    const double *cpoint = geom->GetCurrentPoint();

    //We only exit the loop once we leave the top volume
    while(!(geom->IsOutside())){

        //Get current material
        TGeoMaterial *cmat = CheckMaterial();

        /////////////// Check if the particle is in vacuum or not and propagate accordingly

        if(cmat->GetDensity() == 0){ //Particle is in vacuum (or air, approximatly)

            //std::cout << "Outside scintilator" << endl;

            //Make step to the next boundary and cross it
            geom->FindNextBoundaryAndStep(); //This updates the current position, so cpoint pointer points to the new position

            //std::cout << part->GetEnergy() << std::endl;

            //std::cout << "Step to cross boundary in vacuum: " << step_vac << std::endl;

            //Compute velocity (natural units - c=1)
            v = (double)(part->GetMomentum())/(part->GetEnergy());

            //Compute time
            //t += geom->GetStep()/v; //GetStep returns the step taken with FindNextBoundaryAndStep()

            //Add new point to the particle track
            //track->AddPoint(cpoint[0], cpoint[1], cpoint[2], t);

        } else { //Particle is in some material with specified density

            //std::cout << "Inside scintilator" << endl;

            //Number of arbitrary steps inside material
            int nsteps = 1;

            //The while condition checks whether the next position (if we make the defined step value in the current particle direction)
            //is in the same node/volume of the current position (we leave the loop when the defined step is greater than the distance to the next boundary)
            while(geom->IsSameLocation(cpoint[0] + stepvalue*d[0], cpoint[1] + stepvalue*d[1], cpoint[2] + stepvalue*d[2])){

                //Set arbitrary step - defined in the simulation (by the constructer)
                geom->SetStep(stepvalue);

                //Execute step defined by SetStep (flag = kFALSE means it is an arbitrary step, not limited by geometrical reasons)
                geom->Step(kFALSE); //This updates the current position, so cpoint pointer points to the new position

                //Get velocity before step and update energy and momentum after energy loss
                std::vector<double> aux = Update_EnergyMomentum(stepvalue, part);

                //aux[0] = velocity   aux[1] = energy loss

                //Add step propagation time to current time
                t += double(stepvalue)/aux[0];


                //////////////////////////// Photons generation /////////////////////////////////

                //Depending on the energy lost by the particle during the step
                //there will be a certain number of emited photons by the material (depends on material properties)

                //Get number of photons from Poisson distribution (usually approximatly 1 photon per 100 eV for common scintillators)
                int nphotons = generator->Generate_Photon_Number(aux[1]*10000); //Energy is converted to MeV

                //Save current position in vector (the current geometry position will now be used for the photons)
                //so we need to save the current position of the massive particle to reset the track once we stop propagating the photons
                vector<double> pos {cpoint[0], cpoint[1], cpoint[2]};

                if(nsteps == 100){

                    std::cout << "Number of photons in step " << nsteps << ": " << nphotons << std::endl;

                    //Propagate all emited photons
                    //for(int i = 0; i < nphotons; i++){

                        //generate photon according to Scintillator spectrum
                        Particle* photon = generator->Generate_Photon();

                        //Create secondary Track associated to photon
                        TVirtualGeoTrack* DaughterTrack = track->AddDaughter(1, photon->GetPDG(), photon); //track id, particle pdg, pointer to particle

                        //Assign initial position to the track (the current position of the massive particle)
                        DaughterTrack->AddPoint(pos[0], pos[1], pos[2], t);

                        //Propagate photon
                        PropagatePhoton(DaughterTrack, t);
                    //}
                }

                //Set the current track as the track of the massive particle to continue propagating
                geom->SetCurrentTrack(track_index);

                //Setting current point and direction back to the ones of the massive particle and finding the state
                geom->InitTrack(pos.data(), d.data());

                ++nsteps; //Increse number of steps taken

            } //End of while loop (the particle is now closer to the next boundary than the defined step)

            //Make step to the next boundary and cross it
            geom->FindNextBoundaryAndStep(stepvalue); //This updates the current position, so cpoint pointer points to the new position

            //Get velocity before step and update energy and momentum after energy loss
            std::vector<double> aux2 = Update_EnergyMomentum(geom->GetStep(), part);

            v = aux2[0];

        } // End of propagation inside current material

        //Compute time
        t += geom->GetStep()/v;

        //Add new point to the particle track
        track->AddPoint(cpoint[0], cpoint[1], cpoint[2], t);

    } //End of propagation of the massive particle
}


//////////////////////////////////// Propagate Photon function ///////////////////////////////

void Tracking::PropagatePhoton(TVirtualGeoTrack* track, double t){

    //Set the given track as the current track
    geom->SetCurrentTrack(track);

    //Get photon associated to track
    Particle* part = new Particle(dynamic_cast<Particle*>(track->GetParticle()));

    //Get photon direction
    std::vector<double> d = part->GetDirection();

    //Setting both initial point and direction and finding the state
    geom->InitTrack(track->GetFirstPoint(), d.data());

    //Get new position of the particle after making the step
    const double *cpoint = geom->GetCurrentPoint();

    //Current particle direction
    const double *cdir = geom->GetCurrentDirection();

    //Auxiliary time to propagate the photon
    double t_aux = t;

    //Velocity variable
    double v = 0;

    //Generate random step according to the probability of absorption of the photon
    //(the distance the photon goes through in the material before being absorbed)
    double absorption_step = generator->Generate_Photon_Step();

    std::cout << "Absorption_step: " << absorption_step << std::endl;

    //Distance traveled by the photon in the material
    double total_dist = 0.;

    //We only exit the loop once the photon leaves the top volume
    //or when the distance travelled by the photon in the material is equal to the absorption distance generated (the photon is absorbed)
    while(!(geom->IsOutside())){

        //Get current material
        TGeoMaterial *cmat = CheckMaterial();

        //Vector for new direction (updated after reflection or refraction of light)
        std::vector<double> ndir(3);

        //Get current direction
        const double *cdir = geom->GetCurrentDirection();
        vector<double> dir {cdir[0], cdir[1], cdir[2]};

        /////////////// Check if the particle is in vacuum or not and propagate accordingly

        if(cmat->GetDensity() == 0)
        { //Photon is in vacuum (or air, approximatly)

            v = 1; //Velocity = c in vacuum
            std::cout << "Outside scintilator" << std::endl;

            //Find distance to the next boundary and set step
            geom->FindNextBoundary();

            //Get normal vector to the next boundary
            double* normal = geom->FindNormal(kTRUE);
            vector<double> n(normal,normal+3);

            //Calculate angle between normal vector and incident vector
            double thetai = tools::Angle_Between_Vectors(dir, n);

            //////////////// Check if the photon is going to be reflected or transmited and propagate accordingly

            if(Is_Reflected(thetai))
            { //Photon is reflected

                //Take step to next boundary but dont cross boundary (second flag is kFALSE)
                geom->Step(kTRUE, kFALSE);

                //Get reflected direction
                ndir = tools::Get_Reflected_Dir(dir, n);

            } else
            { //Photon is transmited (refraction according to Snells Law)

                //Take step to next boundary and cross boundary (second flag is kTRUE)
                geom->Step(kTRUE, kTRUE);

                //Get refracted direction
                ndir = tools::Get_Refracted_Dir(dir, n, thetai, 1, 1.58);
            }

            //Set new direction
            geom->SetCurrentDirection(ndir.data());


        }else
        { //Photon is in some material with specified density

            std::cout << "Inside scintilator" << std::endl;

            v = 1/1.58; //Velocity = c/n in some material, where n is the refractive index

            //Find distance to the next boundary and set step
            geom->FindNextBoundary();

            //Get the step taken
            double snext = geom->GetStep();

            std::cout << "Step: " << snext << std::endl;

            //If the absorption distance generated (minus the already travelled distance in this material) is shorter than the distance to the next boundary,
            //the photon is propagated until the point of absorption
            if((absorption_step - total_dist) < snext){

                std::cout << "Absorbed" << std::endl;

                //Set arbitrary step
                geom->SetStep(absorption_step-total_dist);

                //Execute step (flag = kFALSE means it is an arbitrary step, not limited by geometrical reasons)
                geom->Step(kFALSE);

                return; //exit the PropagatePhoton function (the photon is absorbed - no more propagation needed)

            } else {  //if photon reaches the next boundary without being absorbed

                //Get normal vector to the next boundary
                double* normal = geom->FindNormal(kTRUE);
                vector<double> n {normal[0], normal[1], normal[2]};

                //Calculate angle between normal vector and incident vector
                double thetai = tools::Angle_Between_Vectors(dir, n);

                if(Is_Reflected(thetai)){ //Photon is reflected

                    std::cout << "Reflected " << std::endl;

                    //Take step to next boundary but dont cross boundary (second flag is kFALSE)
                    geom->Step(kTRUE, kFALSE);

                    //Get reflected direction
                    ndir = tools::Get_Reflected_Dir(dir, n);

                } else { //Photon is transmited (refraction according to Snells Law)

                    std::cout << "Transmited " << std::endl;

                    //Take step to next boundary and cross boundary (second flag is kTRUE)
                    geom->Step(kTRUE, kTRUE);

                    //Get refracted direction
                    ndir = tools::Get_Refracted_Dir(dir, n, thetai, 1, 1.58);

                }

                //Set new direction
                geom->SetCurrentDirection(ndir.data());

                //Increase distance travelled by the photon in current material
                total_dist += geom->GetStep();
            }

        }

        //Compute time
        t_aux += geom->GetStep()/v;

        //Add new point to the particle track
        track->AddPoint(cpoint[0], cpoint[1], cpoint[2], t_aux);
    }
}


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


////////////////////////////////// Update energy and momentum and return velocity before update //////////////////////////

std::vector<double> Tracking::Update_EnergyMomentum(double step, Particle* part){

    std::vector<double> aux(2); // aux[0] = velocity   aux[1] = energy loss

    //Get energy before step
    double E = part->GetEnergy();

    //Compute velocity before step (natural units - c=1)
    aux[0] = (double)(part->GetMomentum())/E;

    //Calculate energy lost to the material
    aux[1] = BetheBloch(aux[0], step);

    //Update particle energy
    part->ChangeEnergy(E-aux[1]);

    //Update particle momentum
    double p = part->CalculateMomentum(E-aux[1]);
    part->ChangeMomentum(p);

    return aux;
}


////////////////////////////////////// Bethe Bloch function to calculate energy loss in material //////////////////////////////////



//Calculate probability of reflection at boundary
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
bool Tracking::Is_Reflected(double thetai){

    double Reff = FresnelLaw(thetai, 1.58, 1);

    if(generator->Uniform(0,1) < Reff) {

        return true;

    } else {

        return false;
    }
}

////////////////////////////////////////////////// Draw geometry and tracks function //////////////////////////////////////

void Tracking::Draw(){

    TCanvas *c1 = new TCanvas("c1","c1");
    //top->Draw("ogle");
    geom->GetTopVolume()->Draw();
    geom->DrawTracks("/*");
    //geom->AnimateTracks(0, 1000, 200, "/* /G /S");
    c1->SaveAs("Simulation.pdf");
}
