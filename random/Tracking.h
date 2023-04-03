#ifndef __Tracking__
#define __Tracking__

#include "Geometry.h"
#include "TVirtualGeoTrack.h"
#include <vector>
#include "Particle.h"
#include "Generators.h"

class Tracking : public Geometry{

public:
	// constructors and destructor
	Tracking(double radius, double height, double distance, double airgap, double althickness, double step, Generator* g);
	~Tracking() = default;

	//Add track associated to a particle to the geometry
	int AddParticle(int const id, vector<double> x, Particle* particle);

	//Propagate track associated to massive particle
	void Propagate(int track_index);

	//Propagate track associated to optic photon
	void PropagatePhoton(TVirtualGeoTrack* track, double t);

	//Check material where the particle is propagating
	TGeoMaterial* CheckMaterial();

	//Update energy and momentum of particle
	std::vector<double> Update_EnergyMomentum(double step, Particle* part);

	//Calculate energy lost by the particle in small step, when interacting with the material
	double BetheBloch(double v, double step);

	//Photon reflection - gets probability of photon refelection
	double FresnelLaw(double thetai, double n1, double n2);

	//Check if light is reflected or transmitted
	bool Is_Reflected(double thetai);

	void Draw();

protected:

	double stepvalue; //standard propagation step value
	Generator* generator;
};

#endif