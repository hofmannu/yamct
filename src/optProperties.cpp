#include "optProperties.h"

optProperties::optProperties()
{

}


void optProperties::set_mua(const float _mua)
{
	if (_mua >= 0)
	{
		mua = _mua;
	}
	else
	{
		printf("Absorption coeff should be bigger or equal 0");
		throw "invalidValue";
	}

	return;
}

// set scattering coefficient
void optProperties::set_mus(const float _mus)
{
	if (_mus >= 0)
	{
		mus = _mus;
	}
	else
	{
		printf("Scattering coff should be bigger or equal 0");
		throw "invalidValue";
	}
	return;
}

void optProperties::set_g(const float _g)
{
	g = _g;
	return;
}

// set refractive index of medium
void optProperties::set_n(const float _n)
{
	n = _n;
	return;
}