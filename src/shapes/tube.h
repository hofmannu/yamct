/*
	class describing and evaluating the geometrical figure
	of a tube

	Author: Urs Hofmann
	Mail: hofmannu@ethz.ch
	Date: 05.11.2020

	ToDo:
		- add basic functionality like volume and surface area

	Changelog:
		- working version
*/

#ifndef TUBE_H
#define TUBE_H

#include "shape.h"
#include <cstdio>
#include <cstdint>
#include <cmath>

class tube : public shape
{
private:
	float startPos[3] = {0, -1, -1}; // starting position of tube
	float stopPos[3] = {2, 1, 1}; // stop position of tube
	float radius[2] = {0.1, 0.2}; // radius of tube (inner and outer)
	// idx0 is inner radius
	// idx1 is outer radius
public:
	// get functions
	float* get_pradius() {return &radius[0];};
	float* get_pstartPos() {return &startPos[0];};
	float* get_pstopPos() {return &stopPos[0];};

	// returns either inner or outer radius
	float get_radius(const uint8_t idx) const {return radius[idx];};
	float get_iradius() const {return radius[0];};
	float get_oradius() const {return radius[1];};
	float get_startPos(const uint8_t iDim) const {return startPos[iDim];};
	float get_stopPos(const uint8_t iDim) const {return stopPos[iDim];};


	// set functions
	void set_iradius(const float _radius); // define inner radius
	void set_oradius(const float _radius); // define outer radius
	void set_startPos(const float _startPos, const uint8_t iDim);
	void set_stopPos(const float _stopPos, const uint8_t iDim);

	// returns true if position is within tube and false if not
	bool isContained(const float* pos);
};

#endif