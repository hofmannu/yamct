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
#include <nlohmann/json.hpp>

#include "sphere.h"
#include "box.h"
#include "tube.h"

using nlohmann::json;
using namespace std;

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

	vector<sphere> spheres;
	vector<box> boxes;
	vector<tube> tubes;

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
	void generate_volume();

	// adds the shapes to our volume representation
	void add_sphere(sphere* iSphere);
	void add_box(box* iBox);
	void add_tube(tube* iTube);


	// add a new shape definition, do not assign to volume yet
	void new_sphere();
	void new_box();
	void new_tube();

	void delete_sphere(const uint8_t iSphere);
	void delete_box(const uint8_t iBox);
	void delete_tube(const uint8_t iTube);

	// get actual geometrical objects back from vector
	sphere get_sphere(const uint8_t iSphere) const {return spheres[iSphere];};
	box get_box(const uint8_t iBox) const {return boxes[iBox];};
	tube get_tube(const uint8_t iTube) const {return tubes[iTube];};

	sphere* get_psphere(const uint8_t iSphere) {return &spheres[iSphere];};
	box* get_pbox(const uint8_t iBox) {return &boxes[iBox];};
	tube* get_ptube(const uint8_t iTube) {return &tubes[iTube];};

	// set geometrical objects at a certain position
	void set_sphere(const uint8_t iSphere, const sphere sphereObj);
	void set_box(const uint8_t iBox, const box boxObj);
	void set_tube(const uint8_t iTube, const tube tubeObj);

	// get number of geometrical elements of each type
	uint8_t get_nbox() const {return boxes.size();};
	uint8_t get_ntube() const {return tubes.size();};
	uint8_t get_nsphere() const {return spheres.size();};

	float get_volume_voxel(); // returns the volume of each voxel [mm3]
	json get_json(); // returns volume configuration
	void read_json(const json jOpt); 

};

#endif
