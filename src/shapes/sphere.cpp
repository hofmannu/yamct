#include "sphere.h"

// set functions
void sphere::set_radius(const float _radius)
{
	radius = _radius;
	return;
}

void sphere::set_center(const float _center, const uint8_t iDim)
{
	center[iDim] = _center;
	return;
}
// is position part of sphere?
bool sphere::isContained(const float* pos)
{
	// calculate distance between position and center
	float dist = 0;
	float delta;
	for (uint8_t iDim = 0; iDim < 3; iDim++)
	{
		delta = pos[iDim] - center[iDim];
		dist = dist + delta * delta;
	}
	// dist = sqrt(dist); // lets not do this, this is super computationally intense

	bool isContained = (dist <= (radius * radius));

	return isContained;
}