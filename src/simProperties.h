/*
	Properties of used Monte Carlo simulation

	Author: Urs Hofmann
	Mail: hofmannu@ethz.ch
	Date: 05.11.2020

	ToDo:

	Changelog:
		- movd default number of photons to lower values

*/

#ifndef SIMPROPERTIES_H
#define SIMPROPERTIES_H

#include <cstdio>
#include <cstdint>
#include <cmath>

using namespace std;

class simProperties
{
private:
	uint64_t nPhotons = 1e6; // wanted number of photon packages simulated
	uint64_t nPhotonsTrue; // true number of simulated photons
	uint64_t nPPerThread = 10; // number of photons simulated in each thread
	uint64_t threadsPerBlock = 512; // threads per block on GPU
	uint64_t nBlocks; // dependent variable

	void calc_threads();
public:
	simProperties();

	void set_nPhotons(const uint64_t _nPhotons);
	void set_nPPerThread(const uint64_t _nPPerThread);

	uint64_t get_nPhotons() const;
	uint64_t get_nPPerThread() const;
	uint64_t get_nPhotonsTrue() const;
	uint64_t get_threadsPerBlock() const;
	uint64_t get_nBlocks() const;

};
#endif
