/*
	simple class container for optical properties
	Author: Urs Hofmann
	Mail: hofmannu@ethz.ch
	Date: 11.10.2020

	Todo:
		- save and load function
		- automatic determination of tissze props

	Changelog:
		- switched default values to mm
*/
#include <cstdio>

#ifndef OPTPROPERTIES_H
#define OPTPROPERTIES_H

class optProperties
{
private:
	float mua = 0.3; // absorption coefficient [1/mm]
	float mus  = 4.2; // scattering coefficient [1/mm]
	float g = 0.8; // anisotropy of medium
	float n = 1.41; // optical index

	// dependent properties (no set function)
	// float albedo; {
	// albedo = mus / (mus + mua);

public:
	optProperties();
	
	// constant get functions
	float get_mua() const {return mua;};
	float get_mus() const {return mus;};
	float get_g() const {return g;};
	float get_n() const {return n;};
	float get_albedo() const {return (mus / (mus + mua));};

	// pointer return functions
	float* get_pmua() {return &mua;};
	float* get_pmus() {return &mus;};
	float* get_pg() {return &g;};
	float* get_pn() {return &n;};
	// float* get_palbedo() {return &albedo;};


	void set_mua(const float _mua);
	void set_mus(const float _mus);
	void set_g(const float _g);
	void set_n(const float _n);
};

#endif