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
	uint8_t get_tType() const {return tType;};

	void set_tType(const uint8_t _tType);

	bool isContained(const float* pos);

};



#endif