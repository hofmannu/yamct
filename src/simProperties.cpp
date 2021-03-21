#include "simProperties.h"

// constant get properties
uint64_t simProperties::get_nPhotons() const {return nPhotons;}
uint64_t simProperties::get_nPhotonsTrue() const {return nPhotonsTrue;}
uint64_t simProperties::get_threadsPerBlock() const {return threadsPerBlock;}
uint64_t simProperties::get_nBlocks() const {return nBlocks;}

simProperties::simProperties()
{
	calc_threads();
}

void simProperties::set_nPhotons(const uint64_t _nPhotons)
{
	nPhotons = _nPhotons;
	calc_threads();
	return;
}

// set id of gpu to use
void simProperties::set_gpuID(const uint8_t _gpuID)
{
	gpuID = _gpuID;
	return;
}

void simProperties::calc_threads()
{
	nBlocks = ceil( (float) nPhotons / ((float) threadsPerBlock));
	nPhotonsTrue = nBlocks * threadsPerBlock;
	return;
}

void simProperties::set_flagKillBound(const bool _flagKillBound)
{
	flagKillBound = _flagKillBound;
	return;
}