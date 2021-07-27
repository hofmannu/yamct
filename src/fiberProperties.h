/*
	describes the properties of a multimode fiber and their illumination cone

	Author: Urs Hofmann
	Mail: hofmannu@ethz.ch
	Date: 04.11.2020
*/

#ifndef FIBERPROPERTIES_H
#define FIBERPROPERTIES_H

#include <cstring>
#include <iostream>
#include <cmath>

using namespace std;

class fiberProperties
{
private:
	float numAp = 0.2; // numerical aperture of fiber
	float dCore = 0.2; // core diameter of fiber [mm]
	float pos[3] = {0, 0, 0}; // position of fiber input [mm]
	float orientation[3] = {1, 0, 0}; // orientation of fiber input (vector)
	float weight = 1; // factor of photon emissing, e.g 2 means twice as many photons as each 1
	string name;
public:
	
	// class constructor and destructor
	fiberProperties();

	// set functions for fiber properties
	void set_numAp(const float _numAp);
	void set_dCore(const float _dCore);
	void set_rCore(const float _rCore);
	void set_weight(const float _weight);
	void set_name(const string _name);
	void set_pos(const uint8_t iDim, const float _pos);
	void set_orientation(const uint8_t iDim, const float _orientation);

	// constant get functions
	float get_numAp() const {return numAp;};
	float get_dCore() const {return dCore;};
	float get_rCore() const {return dCore / 2;};
	float get_weight() const {return weight;};
	string get_name() const {return name;};
	float get_pos(const uint8_t iDim) const {return pos[iDim];};
	float get_orientation(const uint8_t iDim) const {return orientation[iDim];};

	// pointer return functions
	float* get_pnumAp() {return &numAp;};
	float* get_pdCore() {return &dCore;};
	float* get_ppos() {return &pos[0];};
	float* get_porientation() {return &orientation[0];};
	float* get_pweight() {return &weight;};

	// get functions for dependent properties
	float get_theta(const float nMedium) const; // one sided opening angle
	float get_rSpot(const float nMedium, const float dist) const;

	// not implemented yet
	void saveToFile();
	void saveToFile(const string fileName);
	void readFromFile();
	void readFromFile(const string fileName);
};

#endif
