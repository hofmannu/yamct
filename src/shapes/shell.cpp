#include "shell.h"

// define inner radius of shell
void shell::set_rInner(const float _rInner)
{
	rInner = _rInner;
	return;
}

// define outer radius of shell
void shell::set_rOuter(const float _rOuter)
{
	rOuter = _rOuter;
	return;
}

// define center position along one dimension
void shell::set_center(const uint8_t iDim, const float _center)
{
	center[iDim] = _center;
	return;
}

// define center of shell by passing a pointer
void shell::set_center(const float* _center)
{
	for (uint8_t iDim = 0; iDim < 3; iDim++)
	{
		center[iDim] = _center[iDim];
	}
}

// return either true or false depending if point is contained in shell
bool shell::isContained(const float* pos)
{

	bool amIContained;
	float deltaPos;
	float sumDeltas = 0;
	for (uint8_t iDim = 0; iDim < 3; iDim++)
	{
		deltaPos = pos[iDim] - center[iDim];
		sumDeltas = sumDeltas + deltaPos; 
	}

	if ((sumDeltas >= (rInner * rInner)) && (sumDeltas <= (rOuter * rOuter)))
		amIContained = 1;
	else
		amIContained = 0;

	return amIContained;
}