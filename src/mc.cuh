/*
	main procedure for our monte carlo simulation

	Author: Urs Hofmann
	Mail: hofmannu@ethz.ch
	Date: 04.11.2020

	ToDo:
		- implement oriented fiber options and multiple fibers
		- switch whole heat array to double

	Changelog:
		- switch from radial symmetry approach to cartesian full sized simulation
		- implement fluence mapping on top of heat mapping
		- scaling of heat and fluence to 1 / mm^3 and 1 / mm^2
		- made separate heat bin at heat[nEleme]
		- moved push to heat to more accurate FMA implementation

*/


#ifndef MC_H
#define MC_H

#include "fiberProperties.h"
#include "optProperties.h"
#include "simProperties.h"
#include "optVolume.h" 
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>
#include <math.h>
#include <vector>
#include <H5Cpp.h>


#include "structArgsIn.h"
#include "../lib/vtkwriter/vtkwriter.h"

using namespace std;

class mc
{
private:
	fiberProperties fiber;
	vector<optProperties> tissues; // defines different tissue types
	simProperties sim;
	optVolume volume; // defines the distribution of our tissue types
 
 	// arrays containing our heat map [ix + iy * nx + iz * nx * ny]
	float* heat; // normal heat map
	float* heat_log; // logarthmic scaled heat map
	float* fluence; // fluence map in linear scale
	float* fluence_log; // fluence map in log scale
	bool isHeatAlloc = 0;
	
	// same arrays but this time on our device
	float* heat_dev; // pointer to memory for heat map on device
	float* fluence_dev; // pointer to fluence memory on device
	bool isHeatDevAlloc = 0;
	
	// minimum and maximum values of heat map in lin and log scale
	float maxVal = 0; // maximum value in heat map
	float minVal = 0; // minimum value in heat map
	float maxValLog = 0; // maximum value in logarithmic heat map
	float minValLog = 0; // minimum value in logarthmic heat map

	// min and max values of fluence map in lin and lof scale.
	float maxFluence = 0; // maximum value of fluence
	float minFluence = 0; // minimum value of fluence
	float maxFluenceLog = 0; // maximum log value of fluence
	float minFluenceLog = 0; // minimum log value of fluence

	// things for plotting of cross sections
	uint32_t idxCross[3] = {0, 0, 0}; // current position of our crossection
	float* slice[3]; // slice along each dimension (heat)
	float* sliceLog[3]; // log slice along each dimension (heat)
	float* sliceFluence[3]; // slice of fluence
	float* sliceFluenceLog[3]; // log slice of fluence
	// float* plot[3]; // 1d array in each dimension for plotting
	// float* plotLog[3]; // 1d array in each dimension for plotting 

	uint32_t dims[3]; // dimension of heat map which we get from volume
	float origin[3]; // origin of our volume
	float res[3]; // resolution of our volume

	void init_vars(); // initialize and allocate variables
	
	// void update_plot(const uint8_t iDim, const uint32_t idx0, const uint32_t idx1);
	// void update_plot_log(const uint8_t iDim, const uint32_t idx0, const uint32_t idx1);
	uint32_t idxCrossPlots[3] = {0, 0, 0};
	float* vecs[3]; // vectors along each dimension
	float* plotHeat[3]; // plot vectors along each dimension for heat
	float* plotHeatLog[3]; 
	float* plotFluence[3];
	float* plotFluenceLog[3];

	// void calc_reflectance(const float nWater, const float nTissue); 
	// calculates the tissues interal reflectance vector
	
	void run_sim(); // runs the actual simulation
	void calcMinMax(); // calculate minimum and maximum value of heat map
	void calcLog(); // convert heat and fluence map into log scale
	void update_slice(const uint8_t iDim, const uint32_t idx);
	void update_slice_log(const uint8_t iDim, const uint32_t idx);
	void update_plots(const int32_t* idxPos);

	bool flagDebug = 1; // give more debugging related output in terminal
	bool flagKillBound = 0; // should we kill photon at boundaries

	double simTime = 0; // tiome which was required for last simulation
	bool isDone = 0; // flag showing if simulation is done

	float photonRatio; // ratio of discovered photons in volume

public:

	mc(); // class constructor
	~mc(); // class destructor

	bool get_isDone() const {return isDone;};

	// all our subclass pointer get functions
	fiberProperties* get_pfiber() {return &fiber;};
	simProperties* get_psim() {return &sim;};
	optVolume* get_pvolume() {return &volume;};
	float* get_pheat() {return heat;};
	float* get_pheat_log() {return heat_log;};
	vector<optProperties>* get_ptissues() {return &tissues;};
	float* get_pvec(const uint8_t iDim) {return vecs[iDim];};

	// get maximum and minimum value of mips and mipslog
	float get_maxVal() const {return maxVal;}; // returns maximum value in heat map
	float get_minVal() const {return minVal;}; // returns minimum value in heat map
	float get_maxValLog() const {return maxValLog;}; // returns max log val
	float get_minValLog() const {return minValLog;}; // returns min log val
	
	float get_maxFluence() const {return maxFluence;};
	float get_minFluence() const {return minFluence;};
	float get_maxFluenceLog() const {return maxFluenceLog;};
	float get_minFluenceLog() const {return minFluenceLog;};

	uint32_t get_dim(const uint8_t iDim) const {return dims[iDim];};
	float get_minPos(const uint8_t iDim) const {return origin[iDim];};
	float get_maxPos(const uint8_t iDim) const {return origin[iDim] + res[iDim] * ((float) dims[iDim]);};
	bool get_flagKillBound() const {return flagKillBound;};
	double get_simTime() const {return simTime;};
	float get_photonRatio() const {return photonRatio;};

	void run();

	// returns a crosssection along radial direction
	float* get_pradVec(const float axPos);
	float* get_pradVec(const uint64_t axIdx);
	float* get_pradVec();

	// returns a crosssection along radial direction
	float* get_paxVec(const float radPos);
	float* get_paxVec(const uint64_t radIdx);
	float* get_paxVec();

	bool* get_pflagKillBound() {return &flagKillBound;};
	// returns slice to given position along given dimension
	float* get_slice(const uint8_t iDim, const float pos, const bool flagLog, const bool flagFluence);
	float* get_plot(const uint8_t iDim, const float* pos, const bool flagLog, const bool flagFluence);
	
	// export functions
	bool exportH5(const string filePath); // export file into h5 format
	bool exportH5(); 
	// returns 1 if successful, otherwise 0
	bool exportVtk(const string filePath); // export full heat map to file
};

#endif
