/*
	class representing a box
	Author: Urs Hofmann
	Mail: hofmannu@ethz.ch
	Date: 19.03.2021
*/

#ifndef BOX_H
#define BOX_H

#include <iostream>
#include <math.h>
#include <cinttypes>
#include "shape.h"

class box : public shape
{
private:
	float cornerA[3] = {0, 0, 0}; 
	float cornerB[3] = {1, 1, 1};
public:
	float* get_pcornerA() {return &cornerA[0];};
	float* get_pcornerB() {return &cornerB[0];};
	bool isContained(const float* pos);

	float get_cornerA(const uint8_t iDim) const {return cornerA[iDim];};
	float get_cornerB(const uint8_t iDim) const {return cornerB[iDim];};

	void set_cornerA(const float _value, const uint8_t iDim);
	void set_cornerB(const float _value, const uint8_t iDim);

};

#endif