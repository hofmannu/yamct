/* 
	constant arguments passed to our wondeful CUDA kernel
	used for kernel arguments (constArgsIn) and optical properties

	Author: Urs Hofmann
	Mail: hofmannu@ethz.ch
	Date: 05.11.2020

	ToDo:

	Changelog:
		- added nElements as precalculated value to constArgsIn

*/

#ifndef CONSTARGSIN_H
#define CONSTARGSIN_H

#include <cinttypes>

struct constArgsIn
{
	// properties of output field
	uint32_t dim[3]; // number of elements in field  direction
	uint32_t nElements; // dim[0] * dim[1] * dim[2]
	float res[3]; // resolution / bin size in each direction
	float origin[3]; // first element in each direction
	float maxPos[3]; // stores the maximum position of volume in each dim
	uint8_t bgMaterial = 0; // material of the background


	// angular reflection stuff
	// float critAngle;
	// uint64_t nR;
	// float* Rdev;
	// float dAlpha;
	bool killFlag; // should photons be killed at bounary
};

/*
	struct containing the optical properties of our tissue types
*/
struct optProps
{
	float mu_a; // absorption coefficient in tissue [1/mm]
	float mu_s; // scattering coefficient in tissue [1/mm]
	float g; // anisotropy coefficient
	float n; // refractive index	

	// calculated properties
	float mu_as; // 
	float g2; // nothing but g * g
	float gx2; // nothing but 2 * g
	float albedo; // albedo operator / value
	float albedo1; // 1 - albedo precalculated
};


// properties of fiber output
struct fiberProps
{
	float pos[3]; // position of ecnter of fiber facette
	float orientation[3]; // spatial orientation of optical axis
	float rCore; // core radius of fiber
	float rCore2; // precalculated rCore * rCore
	float na; // numerical aperture of fiber
	float phiMax;
	float randMax; // range of random numbers leading to usage of this fiber
};

#endif
