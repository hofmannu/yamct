/*
	class defining the optical volume we have
	each individual voxel has a byte value which corresponds to the index of the defined
	optical material we are walking through

	Author: Urs Hofmann
	Mail: hofmannu@ethz.ch
	Date: 05.11.2020

	ToDo:

	Changelog:
		- changed default values of origin and resolution

*/

#ifndef OPTVOLUME_H
#define OPTVOLUME_H

#include <cinttypes>
#include <cstdio>
#include <math.h>

#include "sphere.h"
#include "box.h"
#include "tube.h"

class optVolume
{
private:
	float res[3] = {0.01, 0.01, 0.01}; // resolution of volume [mm]
	float lowerCorner[3] = {-0.005, -1.005, -1.005}; // upper corner of volume [mm]
	float upperCorner[3] = {2.005, 1.005, 1.005}; // lower corner of volume [mm]
	uint8_t bgMaterialId = 0; // id of background material
	uint8_t* volumeData; // pointer 
	float* vec[3];
	float* crossSection[3];
	bool isVolumeDataAlloc = 0;

	uint8_t maxMaterial = 0; // highest index of material found in matrix
	uint8_t minMaterial = 0; // minimum index of material found in matrix

	uint32_t dims[3] = {0, 0, 0}; // dimensions of dataset after allocation
	uint32_t nElem = 0; // number of elements in matrix

	float cross[3];

	void calc_dims(); // calculates dimensions and checks validity
public:

	void update_crossSection(const uint8_t iDim, const float crossVal);
	optVolume(); // class constructor
	~optVolume(); // class destructor

	// normal const get functions
	uint8_t get_bgMaterialId() const {return bgMaterialId;};

	uint32_t get_dim(const uint8_t iDim) const {return dims[iDim];};
	float get_res(const uint8_t iDim) const {return res[iDim];};
	uint32_t get_nElem() const {return nElem;};
	float get_min(const uint8_t iDim) const {return lowerCorner[iDim];};
	float get_max(const uint8_t iDim);
	float get_cross(const uint8_t iDim) const {return cross[iDim];};
	uint8_t get_maxMaterial() const {return maxMaterial;};
	uint8_t get_minMaterial() const {return minMaterial;};
	float get_range(const uint8_t iDim) const {return (float) dims[iDim] * res[iDim];};

	// pointer get functions
	float* get_pcrossSection(const uint8_t iDim) {return crossSection[iDim];};
	float* get_pres() {return res;};
	float* get_pupperCorner() {return upperCorner;};
	float* get_plowerCorner() {return lowerCorner;};
	uint8_t* get_pvolumeData() {return volumeData;};
	uint8_t* get_pbgMaterialId() {return &bgMaterialId;};

	// normal set functions
	void set_bgMaterialId(const uint8_t _bgMaterialId);
	void set_cross(const uint8_t iDim, const float crossVal);

	void alloc(); // allocates memory and fill up with background material

	void addSphere(sphere* iSphere);
	void addBox(box* iBox);
	void addTube(tube* iTube);

	float get_volume_voxel(); // returns the volume of each voxel [mm3]


};

#endif
