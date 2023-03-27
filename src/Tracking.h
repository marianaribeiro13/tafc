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

	//Propagate track associated to particle
	void Propagate(int track_index);

	void PropagatePhoton(TVirtualGeoTrack* track, double t);

	//Check material where the particle is propagating
	TGeoMaterial* CheckMaterial();

	//Calculate energy lost by the particle in small step, when interacting with the material
	double BetheBloch(double v, double step);

	//Photon reflection
	double FresnelLaw(double thetai, double n1, double n2);

	//Check if light is reflected or transmitted
	bool Check_Reflection(std::vector<double>& di);


	void Draw();

protected:

	double stepvalue; //standard propagation step value
	Generator* generator;
};

#endif