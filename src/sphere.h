#ifndef SPHERE_H
#define SPHERE_H

#include <math.h>
#include <cinttypes>
#include "shape.h"

class sphere : public shape
{

private:
	float radius; // radius of sphere
	float center[3] = {0.0, 0.0, 0.0}; // center point of sphere
public:

	// get functions
	float* get_pradius() {return &radius;};
	float* get_pcenter() {return &center[0];};
	
	float get_radius() const {return radius;};
	float get_center(const uint8_t iDim) const {return center[iDim];};

	// define radius
	void set_radius(const float _radius); 
	void set_center(const float _center, const uint8_t iDim);
	

	bool isContained(const float* pos);
};

#endif