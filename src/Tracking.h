#ifndef __Tracking__
#define __Tracking__

#include "Geometry.h"

class Tracking : public Geometry{

public:
  // constructors and destructors
  Tracking(Double_t *point, Double_t *dir);
  ~Tracking();
  const Double_t* GetCurPoint();
  const Double_t* GetCurDir();
  TGeoNode* GetCurState();
  void AddTrack(Int_t id, Int_t pdg, TObject* ptrParticle);
  
  //void InitTrack(Double_t *point[3],Double_t *dir[3]);
	/*// getters
	int Ndim() const {return ndim;}
	double T() const {return x[0];}
	double X(int i) const {return x[i+1];}
	double* GetArray() {return x;}

	// operators
	ODEpoint operator*(double) const;
	ODEpoint operator+(const ODEpoint&) const;
	ODEpoint operator-(const ODEpoint&) const;

	void operator=(const ODEpoint&);

	const double& operator[] (int) const;
	double& operator[] (int);

	// print
	friend ostream& operator<<(ostream&, const ODEpoint&);*/

protected:
    
};

#endif
