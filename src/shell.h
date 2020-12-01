/*
	defines the geometry of a shell and maps it to our volume

	Author: Urs Hofmann
	Mail: hofmannu@ethz.ch
	Date: 05.11.2020
*/

#ifndef SHELL_H
#define SHELL_H

#include "shape.h"
#include <cstdio>
#include <cstdint>

class shell : public shape
{

private:
	float center[3]; // center position of shell [mm]
	float rInner = 0; // inner radius of shell [mm]
	float rOuter = 1; // outer radius of shell [mm]

public:
	
	// constant get functions to get values
	float get_rInner() const {return rInner;};
	float get_rOuter() const {return rOuter;};
	float get_center(const uint8_t iDim) const {return center[iDim];};

	// get functions to get pointers to values
	float* get_prInner() {return &rInner;};
	float* get_prOuter() {return &rOuter;};
	float* get_pcenter() {return &center[0];};

	// set functions with control
	void set_rInner(const float _rInner);
	void set_rOuter(const float _rOuter);
	void set_center(const uint8_t iDim, const float _center);
	void set_center(const float* _center);

	bool isContained(const float* pos);

};

#endif