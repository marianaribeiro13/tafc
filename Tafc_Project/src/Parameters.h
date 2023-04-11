#ifndef __Parameters__
#define __Parameters__

struct Parameters{

    ////////////////////// PROPAGATION PARAMETERS ////////////////////////

    int Nmuons_accepted = 2; // Number of muons to propagate (they have to cross both detectors)
    const int N_threads = 2; // Number of threads to use in propagation
    const int muons_per_thread = Nmuons_accepted/N_threads; // Number of accepted muons to propagate on each thread
    
    double step = 0.001; // Tracker step (cm)


    //////////////// TELESCOPE GEOMETRY PARAMETERS (all distances in centimeters) ////////////////////////
    
    const double Radius = 5.0;  //Scintillator radius
	const double Height = 1.0;  //Scintillator height
	const double Distance = 25.;  //Distance between the two scintillators
	const double Thickness = 0.0016; //Aluminium foil thickness
	const double Airgap = .1; //Air gap thickness between the scintillator and the aluminium foil
    
    const double innerradius = Radius + Airgap;
    const double outerradius = Radius + Airgap + Thickness;
	
    int n_SIPM = 4; //Number of SIPMS per scintillator (not being used now)
	const double SIPM_size = .6; //SIPM_size (the detector area is SIPM_size*SIPM_size)
    //Phi angle of each SIPM placed on center of lateral scintillator wall (size of the vector is the number of SIPMs per scintillator)
    std::vector <double> SIPM_angles {0., M_PI/2, M_PI, 3*M_PI/2}; 
	
    const double SIPM_angle = 2*M_PI / n_SIPM;
	const double SIPM_alpha = 0.5*SIPM_size/(SIPM_angle*Radius);
    const double SIPM_phi_range = SIPM_size/Radius; //Approximate SIPM phi range (assumes SIPM size much smaller than scintillator radius)
    

};

#endif