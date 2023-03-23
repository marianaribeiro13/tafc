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
	Tracking(double distance, double step, Generator* g);
	~Tracking() = default;

	//Add track associated to a particle to the geometry
	int AddParticle(int const id, vector<double> x, Particle* particle);

	//Propagate track associated to particle
	void Propagate(int track_index);

	//Check material where the particle is propagating
	TGeoMaterial* CheckMaterial();

	//Calculate energy lost by the particle in small step, when interacting with the material
	double BetheBloch(double v, double step);

	Particle* GenerateCosmicMuon();

	void Draw();

protected:

	double stepvalue; //standard propagation step value
	Generator* generator;
};

#endif