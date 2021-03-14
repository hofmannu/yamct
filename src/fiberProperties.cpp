#include "fiberProperties.h"

fiberProperties::fiberProperties()
{
		name = "default";

}

void fiberProperties::set_numAp(const float _numAp)
{
	if (_numAp >= 0)
		numAp = _numAp;
	else
	{
		printf("Numerical aperture must be equal or bigger 0\n");
		throw "invalidValue";
	}
	return;
}

void fiberProperties::set_dCore(const float _dCore)
{
	if (dCore >= 0)
	{
		dCore = _dCore;
	}
	else
	{
		printf("Core diameter must be bigger then 0\n");
		throw "invalidValue";
	}
	return;
}

void fiberProperties::set_rCore(const float _rCore)
{
	dCore = 2 * _rCore;
	return;
}

void fiberProperties::set_weight(const float _weight)
{
	weight = _weight;
	return;
}

void fiberProperties::set_name(const string _name)
{
	name = _name;
	return;
}

float fiberProperties::get_theta(const float nMedium) const
{
	float theta = asin(numAp / nMedium);
	return theta;
}

float fiberProperties::get_rSpot(const float nMedium, const float dist) const
{
	float rSpot = get_rCore() + tan(get_theta(nMedium)) * dist;
	return rSpot;
}

void fiberProperties::set_pos(const uint8_t iDim, const float _pos)
{
	pos[iDim] = _pos;
	return;
}

void fiberProperties::set_orientation(const uint8_t iDim, const float _orientation)
{
	orientation[iDim] = _orientation;
	return;
}