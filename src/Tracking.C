#include "Tracking.h"
#include "TCanvas.h"
#include "TMath.h"


////////////////////////// Constructor ////////////////////////////////////////
Tracking::Tracking(double distance, double step, Generator* g) : Geometry(distance), stepvalue(step){
    generator = g;
}

//Adds track associated to particle (give initial position)
int Tracking::AddParticle(int const id, vector<double> x, Particle* particle){

    //Create Track associated to particle
    int track_index = geom->AddTrack(id, particle->GetPDG(), particle); //track id, particle pdg, pointer to particle

    //Get created track
    TVirtualGeoTrack* track = geom->GetTrack(track_index);

    //Assign initial position to the track
    track->AddPoint(x[0], x[1], x[2], 0);

    return track_index;
}


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
            double snext  = gGeoManager->GetStep();

            //Compute velocity
            double v = (double)(part->GetMomentum())/(part->GetEnergy());

            //Compute time
            t += double(snext)/v;

            //Get new position of the particle after making the step
            const double *cpoint = geom->GetCurrentPoint();

            //Add new point to the particle track
            track->AddPoint(*cpoint, *(cpoint+1), *(cpoint+2), t);

        } else {

            std::cout << "Inside scintilator" << endl;

            //Get current particle position in the geometry
            const double *ipoint = geom->GetCurrentPoint();
            //Get current particle direction of propagation
            const double *idir = geom->GetCurrentDirection();

            int nsteps = 0;

            while(geom->IsSameLocation(*ipoint + nsteps*stepvalue*(*idir), *(ipoint+1) + nsteps*stepvalue*(*(idir+1)), *(ipoint+2) + nsteps*stepvalue*(*(idir+2)))){
                
                //std::cout << "Step " << nsteps << std::endl;
                //Set arbitrary step
                geom->SetStep(stepvalue);

                //Execute step (flag = kFALSE means it is an arbitrary step, not limited by geometrical reasons)
                geom->Step(kFALSE);

                //Compute particle velocity
                double v = (double)(part->GetMomentum())/(part->GetEnergy());

                //std::cout << v << std::endl;

                //Compute propagation time
                t += double(stepvalue)/v;

                //Get new position of the particle after making the step
                const double *npoint = geom->GetCurrentPoint();

                //Add new point to the particle track
                track->AddPoint(*npoint, *(npoint+1), *(npoint+2),t);

                //Get energy before step
                double E = part->GetEnergy();

                std::cout << "Energy step " << nsteps << ": " << E << std::endl;

                //Calculate energy lost to the material
                double dE = BetheBloch(v, stepvalue);

                //Update particle energy
                part->ChangeEnergy(E-dE);

                //std::cout << part->GetEnergy() << std::endl;

                //Update particle momentum
                double p = part->CalculateMomentum(E-dE);
                part->ChangeMomentum(p);

                nsteps++;
            }
                
            //Make step to the next boundary and cross it
            geom->FindNextBoundaryAndStep();

            //Get the step taken
            double snext  = gGeoManager->GetStep();

            //Compute velocity
            double v = (double)(part->GetMomentum())/(part->GetEnergy());

            //Compute time
            t += double(snext)/v;

            //Get new position of the particle after making the step
            const double *cpoint = geom->GetCurrentPoint();

            //Add new point to the particle track
            track->AddPoint(*cpoint, *(cpoint+1), *(cpoint+2), t);

            //Get energy before step
            double E = part->GetEnergy();

            //Calculate energy lost to the material
            double dE = BetheBloch(v, snext);

            std::cout << "Step to cross boundary: " << snext << std::endl;

            //Update particle energy
            part->ChangeEnergy(E-dE);

            //Update particle momentum
            double p = part->CalculateMomentum(E-dE);
            part->ChangeMomentum(p);

        }
    }
}

TGeoMaterial* Tracking::CheckMaterial(){

    //Get current node
    TGeoNode *cnode = geom->GetCurrentNode();

    //Get current volume
    TGeoVolume *cvol = cnode->GetVolume();

    //Get current material
    TGeoMaterial *cmat = cvol->GetMedium()->GetMaterial();

    return cmat;
}

double Tracking::BetheBloch(double v, double step){

    //NATURAL
    //Double_t alpha=0.0072974;
    //Double_t I1=64.70713358e-6; //MeV
    //Double_t n=2.5587248e18; //MeV
    Int_t Z=-1;
    //Double_t mass_muon=105.6583755; //MeV
    //Double_t mass_electron=0.51099895; //MeV
    //Double_t dE=0;
    Double_t dx=1.e-4;

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

    double dE_SI = dEdx_SI * step;

    return dE_SI;
}

Particle* Tracking::GenerateCosmicMuon(){

    //2D Spectrum of the muon (x[0] is the momentum and x[1] is the incident angle)
   /* auto f = [](double *x,double *par)
    {
        return pow(cos(x[1]),3)*0.00253*pow(x[0]*cos(x[1]),-(0.2455+1.288*log10(x[0]*cos(x[1]))-0.255*pow(log10(x[0]*cos(x[1])),2)+0.0209*pow(log10(x[0]*cos(x[1])),3) ));
    };

    //Create TF1 function with the 2D spectrum of the muon
    TF1 *F = new TF1("f",f);*/

    //Generate random momentum and incident angle according to 2D pdf using acceptance-rejection
    std::vector<double> aux = generator->Random_Distribution_2D(generator->GetMomentumDistribution(),1,2000,0, TMath::Pi()/2,generator->GetMomentumDistribution()->GetMaximum());

    //Get momentum in MeV
    //double Momentum = aux[0]*1000;

    //Get incident angle
    //double theta = aux[1];

    //Get muon direction from the incident angle and a random generation on xy plane
    //vector<double> d = generator->Generate_Direction_From_Theta(aux[1]);

    //Create muon (pdg = 13)
    Particle* muon = new Particle(13, aux[0]*1000, generator->Generate_Direction_From_Theta(aux[1]));

    
    //std::cout << muon->GetEnergy() << std::endl;
    //std::cout << "Momentum: " << aux[0]*1000 << std::endl;

    return muon;
}

void Tracking::Draw(){

    TCanvas *c1 = new TCanvas("c1","c1");
    c1->cd();
    //top->Draw("ogle");
    top->Draw();
    geom->DrawTracks("/*");
    //geom->AnimateTracks(0, 1000, 200, "/* /G /S");
    c1->SaveAs("Drawing.pdf");
}
