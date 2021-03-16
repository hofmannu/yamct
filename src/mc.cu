#include "mc.cuh"
#include "quaternion_toolbox.cu"

// scatters the photon into a new direction depending on currently employed 

// perfrom quaternion rotation
//   - r is normalized vector of rotation
//   - alpha is angle
//   - xyz is point to rotate, will be modified directly

// returns the optical properties for a given position
__device__ __inline__ static optProps getOptProps(
	vec3& pos, // requested position as 3 element vector
	const constArgsIn& inArgs, // constant input arguments to kernel
	const uint8_t* tissueArray, // array defining tissue properties
	const optProps* tissues) // optical properties of tissue
{
	// get position in grid
	int32_t posIdx [3];
	uint8_t tissueId = inArgs.bgMaterial;	

	posIdx[0] = (pos.x - inArgs.origin[0]) / inArgs.res[0];
	

	// if we are outside tissue range, simply assign bgMaterial
	if ((posIdx[0] >= 0) && (posIdx[0] < inArgs.dim[0]))
	{
		posIdx[1] = (pos.y - inArgs.origin[1]) / inArgs.res[1];
		if ((posIdx[1] >= 0) && (posIdx[1] < inArgs.dim[1]))
		{
			posIdx[2] = (pos.z - inArgs.origin[2]) / inArgs.res[2];
			if ((posIdx[2] >= 0) && (posIdx[2] < inArgs.dim[2]))
			{
				const uint32_t iTissue = posIdx[0] + inArgs.dim[0] * 
					(posIdx[1] + inArgs.dim[1] * posIdx[2]);
				// get index of tissue
				tissueId = tissueArray[iTissue];
			}
		}
	}

	// return corresponding tissue
	return tissues[tissueId];
}

// returns random angle required for scattering of photon into a new direction
__device__  __inline__ static float getsigma(
	const float& g, // anisotropy coefficient
	const float& gx2, // anisotropy coefficient times two
	const float& g2, // anisotropy coefficient power two
	curandState& state)
{

	float sigma;
	// scatter photon
	if (g == 0)
	{
		// sigma = acosf(__fmul_rn(2.0, curand_uniform(&state)) - 1.0); 
		sigma = acosf(fmaf(2.0, curand_uniform(&state), -1.0));
		// cos(sigma) = 2 * RND - 1
	}else{
		// g - anisotroy coefficient, g2 = g^2, gx2 = 2 * g
		const float temp = __fdividef(
			1 - g2, 
			fmaf(gx2, curand_uniform(&state), 1 - g)
		);
		// sigma = acosf(__fdividef((1 + g2 - __fmul_rn(temp, temp)), gx2));
		sigma = acosf(__fdividef(fmaf(-temp, temp, 1 + g2), gx2));
	}

	return sigma;
}

// normalize our velocity vector if it really got out of boounds
__device__  __inline__ static void normalize(vec3* vector)
{
	// const float normVal = rnorm3df(vector->x, vector->y, vector->z);
	float normVal = fmaf(vector->z, vector->z, fmaf(vector->y, vector->y, vector->x * vector->x));

	if ((normVal < 0.95) || (normVal > 1.05))
	{
		normVal = sqrtf(normVal);
		printf("velocity is off by a lot: %f!\n", normVal);
		if (normVal == 0)
		{
			printf("Norm 0 is a crime, you rude human!\n");
		}
		else
		{
			vector->x /= normVal;
			vector->y /= normVal;
			vector->z /= normVal;
		}
	}

	return;
}

// implemented fmaf for increased accuracy
__device__  __inline__ static void scatter(
	vec3* vel, // most recent directory
	curandState& state, // random number generator
	const optProps tissueProps) // properties of tissue we are currently in
{
	// generate rotation angles for scattering
	const float sigma = getsigma(
		tissueProps.g, tissueProps.gx2, tissueProps.g2, state);
	const float phi = __fmul_rn(__fmul_rn(2.0, M_PI), curand_uniform(&state));
	vec3 b; // rotation vector oriented perpendicular to vel and new direction

	if ((fabs(vel->x) >= fabs(vel->y)) && (fabs(vel->x) >= fabs(vel->z))){
		// photon is mostly moving in x, so rectangular vector lies in y z plane

		// b.y = abs(vx / sqrt(1 + vx^2 + 2 * vy * vz))
		b.y = fabs(__fdividef(vel->x, 
				__fsqrt_rn(fmaf(2.0, __fmul_rn(vel->y, vel->z), fmaf(vel->x, vel->x, 1)))
			)
		);
		b.z = b.y;
		// b.x = -b.y / vel.x * (vel.y + vel.z);
		b.x = -__fmul_rn( __fdividef(b.y, vel->x), (vel->y + vel->z));
		//printf("A - Initial deflection vector b: %f, %f, %f\n", b.x, b.y, b.z);
	}
	else if ((fabs(vel->y) >= fabs(vel->x)) && (fabs(vel->y) >= fabs(vel->z))){
		// photon is mostly moving in y, so rectangular vector lies in x y plane
		// b.x = abs(vy / sqrt(1 + vy^2 + 2 * vx * vz))
		b.x = fabs(__fdividef(vel->y,
				__fsqrt_rn(fmaf(2.0, __fmul_rn(vel->x, vel->z), fmaf(vel->y, vel->y, 1)))
			)
		);
		b.z = b.x;
		// b.y = -b.x / vy * (vx + vz)
		b.y = -__fmul_rn(__fdividef(b.x, vel->y), (vel->x + vel->z));
		//printf("B - Initial deflection vector b: %f, %f, %f\n", b.x, b.y, b.z);
	}
	else if ((fabs(vel->z) >= fabs(vel->x)) && (fabs(vel->z) >= fabs(vel->y))){
		// photon is mostly moving in z, so rectangular vector lies in x y plane
		// b.x = abs(vz / sqrt(1 + vz^2 + 2 * vx * vy))
		b.x = fabs(__fdividef(vel->z,
				__fsqrt_rn(fmaf(2.0, __fmul_rn(vel->x, vel->y), fmaf(vel->z, vel->z, 1)))
			)
		);
		b.y = b.x;
		b.z = -__fmul_rn(__fdividef(b.x, vel->z), (vel->x + vel->y));
		//printf("C - Initial deflection vector b: %f, %f, %f\n", b.x, b.y, b.z);
	}
	else
	{
		printf("something here does not work at all\n");
		printf("vel vector = [%f, %f, %f]\n", vel->x, vel->y, vel->z);
	}


	// printf("rot vec length: %f, vel before: %f\n", getLengthVec(b), getLengthVec(vel));
	quaternion_rot(vel, phi, &b); // rotate b around vel to generate deflection vector	
	quaternion_rot(&b, sigma, vel); // rotate vel around rotated b by sigma
	normalize(vel);
	// printf("rot vec length: %f, vel after: %f\n", getLengthVec(b), getLengthVec(vel));

	return;
}

/*
	launch a new photon from the fiber tip into our space
	updated to 3D version already
	launch(weight, pos, vel, states[tid], fiber, tissueTypes);

	implemented fmaf for icnreased accuracy
*/
__device__  __inline__ static void launch(
	vec3& pos, // position of photon in [x, y, z] in [m]
	vec3* vel, // directory of movement in [x, y, z]
	curandState& state, // random number generator
	const fiberProps* fibers)
{

	// check which fiber we are using
	fiberProps currFiber; // temporary variable for fiber properties
	const float randFiber = curand_uniform(&state); // random number
	uint8_t idx = 0;
	do
	{

		currFiber = fibers[idx];
		idx++;
	}while (currFiber.randMax < randFiber);

	// generate random start positions in y z plane until in rCore
	vec3 fiberOff;
	do{
		// ((rand * 2.0) - 1) * rCore
		fiberOff.z = fmaf(curand_uniform(&state), 2.0, -1) * currFiber.rCore; 
		fiberOff.y = fmaf(curand_uniform(&state), 2.0, - 1) * currFiber.rCore;
	}while(currFiber.rCore2 < (fmaf(fiberOff.z, fiberOff.z, fiberOff.y * fiberOff.y)));
	fiberOff.x = 0;

	float dirY, dirZ, rUnit;
	do{
	  dirY = fmaf(curand_uniform(&state), 2.0, - 1);
	  dirZ = fmaf(curand_uniform(&state), 2.0, - 1);
    rUnit = fmaf(dirY, dirY, dirZ * dirZ);
  }while(rUnit > 1);
   
  const float sinPhiMax = sinf(currFiber.phiMax);
  vel->y = dirY * sinPhiMax;
  vel->z = dirZ * sinPhiMax;
    
  rUnit = fmaf(vel->y, vel->y, vel->z * vel->z); // rUint = y * y + z * z
  vel->x = sqrtf(1 - rUnit);

	
  // construct rotation vector for quaternion with length 1
  // bug here!!
  if (!((currFiber.orientation[1] == 0) && (currFiber.orientation[2] == 0)))
  {
	  vec3 rotVector;
	  const float normVec = sqrtf(
	  	fmaf(currFiber.orientation[1], currFiber.orientation[1], 
	  	currFiber.orientation[2] * currFiber.orientation[2]));
	  rotVector.x = 0;
	  rotVector.y = -currFiber.orientation[2] / normVec;
	  rotVector.z = currFiber.orientation[1] / normVec;

	  printf("Initial vel vector: %f, %f, %f\n", vel->x, vel->y, vel->z);
	  // calulate rotation angle about the defined rotation vector
	  const float sigma = acos(currFiber.orientation[0]);
		quaternion_rot(&rotVector, sigma, vel); // rotate direction vector
	  quaternion_rot(&rotVector, sigma, &fiberOff); // rotate facette
	}

	// update position to fiber t
	pos.x = currFiber.pos[0] + fiberOff.x;
	pos.y = currFiber.pos[1] + fiberOff.y;
	pos.z = currFiber.pos[2] + fiberOff.z;

	return;
}

// return the minimum of two values
__device__  __inline__ static float min2(const float a, const float b) 
{
	const float m = (a >= b) ? a : b;
	return m;
}

// return the minnimum of two values
__device__  __inline__ static float min3(const float a, const float b, const float c)
{
	float m;
	if (a <=  min2(b, c))
	    m = a;
	else if (b <= min2(a, c))
	    m = b;
	else
	    m = c;
	return m;
}

// check if two positions are contained within the same voxel
__device__ __inline__ static bool sameVoxel(
	const vec3& posA, 
	const vec3& posB, 
	const constArgsIn& inArgs)
{	
	bool sv = 0; // initially we are pessimistic (not same voxel)
  const int32_t idxAX = posA.x / inArgs.res[0];
  const int32_t idxBX = posB.x / inArgs.res[0];
  if (idxAX == idxBX)
  {
  	const int32_t idxAY = posA.y / inArgs.res[1];
  	const int32_t idxBY = posB.y / inArgs.res[1];
  	if (idxAY == idxBY)
  	{
  		const int32_t idxAZ = posA.z / inArgs.res[2];
  		const int32_t idxBZ = posB.z / inArgs.res[2];
  		if (idxAZ == idxBZ)
  		{
  			sv = 1; // only if all three inidices are the same we continue
  		}
  	}
  }
  return sv;
}

// Implementation of FMA as atomic operation for floats
// z = z + x * y
__device__ __inline__ float atomicFMA(float* address, const float x, const float y)
{
  int *address_as_ull = (int*) address;
	int old = *address_as_ull;
	int assumed;

	do{
		assumed = old;
	  old = atomicCAS(
	  	address_as_ull, 
	  	assumed,
	  	__float_as_int(fmaf(x, y, __int_as_float(assumed))));

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
  } while (assumed != old);

  return __float_as_int(old);
}

// push value over to heat matrix
__device__ __inline__ void pushToHeatFMA(
	const vec3& pos, // position of where our absorption event is happening 
	const float weight, // amount of absorption event
	const float scaleH, // weight of the photon stored in heat map
	float* heat, // matrix containing the heat
	const constArgsIn& inArgs) // constant input arguments for kernel
{
	// convert position into index
	bool inside = 0;

	// calculate index in x, y, z direction
	const int64_t ix = (pos.x - inArgs.origin[0]) / inArgs.res[0];
	if ((ix >= 0) && (ix < inArgs.dim[0]))
	{
		const int64_t iy = (pos.y - inArgs.origin[1]) / inArgs.res[1];
		if ((iy >= 0) && (iy < inArgs.dim[1]))
		{
			const int64_t iz = (pos.z - inArgs.origin[2]) / inArgs.res[2];
			if ((iz >= 0) && (iz < inArgs.dim[2]))
			{
				// assign absorption event to heat map
				int64_t idx = ix + inArgs.dim[0] * (iy + iz * inArgs.dim[1]);
				atomicFMA(&heat[idx], weight, scaleH);
				inside = 1;
			}
		}
	}

	if (inside == 0) // if we are outside, add to container bin instead
	 	atomicFMA(&heat[inArgs.dim[0] * inArgs.dim[1] * inArgs.dim[2]], weight, scaleH);

	return;
}


/* 
	find time required to reach closest voxel face in our direction of movement
*/
__device__ __inline__ float findDistFace(
	const vec3& posA, 
	const constArgsIn& inArgs,
	const vec3& vel)
{	
	// find lower side of voxel in index offset w/ origin  = 0
  const int ix1 = floor(posA.x / inArgs.res[0]);
  const int iy1 = floor(posA.y / inArgs.res[1]);
  const int iz1 = floor(posA.z / inArgs.res[2]);
  
  // depending on velocity direction define if upper or lower side meets path 
  const int ix2 = (vel.x >= 0) ? (ix1 + 1) : (ix1);
  const int iy2 = (vel.y >= 0) ? (iy1 + 1) : (iy1);
  const int iz2 = (vel.z >= 0) ? (iz1 + 1) : (iz1);
  
  // convert closest face into distance
  const float deltaX = fmaf((float) ix2, inArgs.res[0], -posA.x);
	const float deltaY = fmaf((float) iy2, inArgs.res[1], -posA.y);
	const float deltaZ = fmaf((float) iz2, inArgs.res[2], -posA.z);		

	// convert through velocity into "time"
	const float tx = fabs(deltaX / vel.x);
  const float ty = fabs(deltaY / vel.y);
  const float tz = fabs(deltaZ / vel.z);
    
  const float tmin = min3(tx, ty, tz);

  return tmin;
}

// move our photons position along vel vector by step
__device__ __inline__ void move(vec3& pos, const vec3& vel, const float step)
{	
	pos.x = fmaf(vel.x, step, pos.x); // pos.x += vel.x * step;
	pos.y = fmaf(vel.y, step, pos.y); // pos.y += vel.y * step;
	pos.z = fmaf(vel.z, step, pos.z); // pos.z += vel.z * step;
	return;
}

// function sets weight to 0 if photon is beyond boundaries and 
// kills it if this is the case
__device__ __inline__ void checkBoundary(
	const vec3& pos, // current position of photon 
	const constArgsIn& inArgs, // constant input arguments
	float& weight, // weight of photon
	float* heat_dev)
{
	if ((pos.x < inArgs.origin[0]) || (pos.x > inArgs.maxPos[0]))
	{
		if ((pos.y < inArgs.origin[1]) || (pos.y > inArgs.maxPos[1]))
		{
			if ((pos.z < inArgs.origin[2]) || (pos.z > inArgs.maxPos[2]))
			{				
				atomicAdd(&heat_dev[inArgs.nElements], weight);
				weight = 0;
			}
		}
	}
	
	return;
}

// cuda kernel definition
__global__ void simPhoton
(
	float * heat_dev, // output matrix to which we write our absorption
	const constArgsIn* inArgsPtr, // constant simulation parameters
	const uint8_t* tissueTypes, // vector containing current tissue props
	const optProps* optProp_dev,
	const fiberProps* fiberProp_dev) // vector containing fiber properties
{
	const constArgsIn inArgs = inArgsPtr[0];

	// temp variables for scattering
	const unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x; 
	const float ls = 1e-3;// minimum travel distance [mm]
		
	curandState cuState; // generate and initialize random number generator 
	curand_init(tid, 0, 0, &cuState); // Initialize CURAND for this thread

	// generate starting position
	vec3 pos; // position of photon [m]
	vec3 vel; // directivity of photon propagation in tissue, always length == 1
	vec3 tempPos;
	optProps currProps; // temporary variable for optical properties of tissue

	float s; // actual distance traveled by photon [mm]
	float sleft; // unitless hop distance
	// to get from unitless to unitfull hop distance: s_uint = s_unitless / mus
	bool sv; // flag defining if we move within the same voxel
		
	float weight = 1.0;
	launch(pos, &vel, cuState, fiberProp_dev);
	currProps = getOptProps(pos, inArgs, tissueTypes, optProp_dev);

	// here starts the intense loop where calculation speed is critical
	unsigned int iScatter = 0;
	while(weight > 0)
	{
		sleft = -__logf(curand_uniform(&cuState)); // dimensionless pathlength
		do
		{
			s = sleft / currProps.mu_s;
			tempPos = pos; // get position estimate oif within same material
			move(tempPos, vel, s); // tempPos = tempPos + vel * s
			
			sv = sameVoxel(pos, tempPos, inArgs); // check if temporary pos are in same voxel
			if (sv) // if our photon is in same voxel
			{ 
				pos = tempPos; // update current position
				const float scaleH = -expm1f(-currProps.mu_a * s); // 1 - exp(-mua * s)
				pushToHeatFMA(pos, weight, scaleH, heat_dev, inArgs);
				weight = fmaf(-weight, scaleH, weight); // weight = weight * (1 - scaleH);
				sleft = 0; // update sleft
			}
			else // photon has crossed voxel boundary
			{
				s = ls + findDistFace(pos, inArgs, vel); // actual travel distance [mm]
				const float scaleH = -expm1f(-currProps.mu_a * s); // 1 - exp(-mua * s)
				pushToHeatFMA(pos, weight, scaleH, heat_dev, inArgs);
				weight = fmaf(-weight, scaleH, weight); // weight = weight * (1 - scaleH);
				sleft = fmaf(-s, currProps.mu_s, sleft); // fmaf(x, y, z) --> x * y + z
				move(pos, vel, s); // update positions

				if (sleft <= 0.01) // check if below littelest unitless step
					sleft = 0;
				
				// update optical properties of tissue
				currProps = getOptProps(pos, inArgs, tissueTypes, optProp_dev);
			}
								
			if (inArgs.killFlag) // if kill flag enabled, check if still in range
				checkBoundary(pos, inArgs, weight, heat_dev);

			if (weight < 1e-3) // play roulette with photon
			{
				if (curand_uniform(&cuState) > 0.1)
					weight = 0; // kill photon
				else
					weight = __fmul_rn(weight, 10);
			}
		
		}while((sleft > 0) && (weight > 0)); // iteratively move and absorb in here
		scatter(&vel, cuState, currProps); // scatter photon if required
		// alterantive implementation as scatterLW (contains steps)
	}
	return;
}


// class destructor mostly to free up all the memory
mc::~mc()
{
	// free volumetric data for our heat map
	if (isHeatAlloc)
	{
		delete[] heat;
		delete[] heat_log;
		delete[] fluence;
		delete[] fluence_log;
		// free slices and plots
		for (uint8_t iDim = 0; iDim < 3; iDim++)
		{
			delete[] slice[iDim];
			delete[] sliceLog[iDim];
			delete[] sliceFluence[iDim];
			delete[] sliceFluenceLog[iDim];
			delete[] vecs[iDim];
			delete[] plotHeat[iDim];
			delete[] plotHeatLog[iDim];
			delete[] plotFluence[iDim];
			delete[] plotFluenceLog[iDim];
		}
	}

	if (isHeatDevAlloc)
	{
		cudaFree(heat_dev);
		cudaFree(fluence_dev);
	}

}

// class constructor
mc::mc()
{
	// generate default material 0 as water
	optProperties water;
	water.set_mua(0.4); // 1/mm
	water.set_mus(2);
	water.set_n(1.31);
	water.set_g(0.9); // does not matter as long as we have no scattering

	tissues.push_back(water);

	fiberProperties startFiber;
	fibers.push_back(startFiber); 
}

float* mc::get_slice(
	const uint8_t iDim, // dimension which we cut 
	const float pos, // position where we cut
	const bool flagLog, // flag for log scale
	const bool flagFluence) // flag for fluence scale
{
	// convert position into index
	int32_t idx = (pos - origin[iDim]) / res[iDim];
	
	if (idx < 0) // if below 0, limit to 0
		idx = 0;

	if (idx >= dims[iDim]) // if beyond upper range, also limit
		idx = dims[iDim] - 1;

	// check if idx is same as last time, then simply return, otherwise
	// update slice
	if (idx != idxCross[iDim])
	{
		update_slice_log(iDim, idx);
		update_slice(iDim, idx);
		idxCross[iDim] = idx;
	}

	float* retSlice;

	if (flagFluence)
	{
		if (flagLog)
		{
			retSlice = sliceFluenceLog[iDim];
		}
		else
		{
			retSlice = sliceFluence[iDim];
		}
	}
	else
	{
		if (flagLog)
		{
			retSlice = sliceLog[iDim];
		}
		else
		{
			retSlice = slice[iDim];
		}
	}

	return retSlice;
}

float* mc::get_plot(
	const uint8_t iDim, // dimension along which we request
	const float* pos, // three element vector describing center position
	const bool flagLog,
	const bool flagFluence)
{
	int32_t idx[3];
	bool flagChange = 0;

	// convert position in all three dimensions into an index
	for (uint8_t iDim = 0; iDim < 3; iDim++)
	{
		idx[iDim] = (pos[iDim] - origin[iDim]) / res[iDim];
		if (idx[iDim] < 0)
			idx[iDim] = 0;

		if (idx[iDim] >= dims[iDim])
			idx[iDim] = dims[iDim] - 1;

		if (idx[iDim] != idxCrossPlots[iDim])
			flagChange = 1;
	}

	// if any position changed we update all plots
	if (flagChange)
	{
		update_plots(&idx[0]);
		for (uint8_t iDim = 0; iDim < 3; iDim++)
		{
			idxCrossPlots[iDim] = idx[iDim];
		}
	}

	float* returnedPlot;

	if (flagFluence) // here we are supposed to return fluence
	{
		if (flagLog) // fluence in log
		{
			returnedPlot = plotFluenceLog[iDim];
		}
		else // fluence in normal
		{
			returnedPlot = plotFluence[iDim];
		}
	}
	else // now we push back with heat
	{
		if (flagLog) // heat in log scale
		{
			returnedPlot = plotHeatLog[iDim];
		}
		else // heat in normal scale
		{
			returnedPlot = plotHeat[iDim];
		}
	}	
	return returnedPlot;
}

// updates all plots along all dimensions
void mc::update_plots(const int32_t* idxPos)
{

	uint32_t idxIn; // iX + nX * (iY * iZ * nY)
	float* currPlot;

	for (uint32_t iX = 0; iX < dims[0]; iX++)
	{
		idxIn = iX + dims[0] * (idxPos[1] + idxPos[2] * dims[1]);
		
		currPlot = plotHeat[0];
		currPlot[iX] = heat[idxIn];

		currPlot = plotHeatLog[0];
		currPlot[iX] = heat_log[idxIn];

		currPlot = plotFluence[0];
		currPlot[iX] = fluence[idxIn];

		currPlot = plotFluenceLog[0];
		currPlot[iX] = fluence_log[idxIn];
	}

	for (uint32_t iY = 0; iY < dims[1]; iY++)
	{
		idxIn = idxPos[0] + dims[0] * (iY + idxPos[2] * dims[1]);

		currPlot = plotHeat[1];
		currPlot[iY] = heat[idxIn];

		currPlot = plotHeatLog[1];
		currPlot[iY] = heat_log[idxIn];

		currPlot = plotFluence[1];
		currPlot[iY] = fluence[idxIn];

		currPlot = plotFluenceLog[1];
		currPlot[iY] = fluence_log[idxIn];
	}

	for (uint32_t iZ = 0; iZ < dims[2]; iZ++)
	{
		idxIn = idxPos[0] + dims[0] * (idxPos[1] + iZ * dims[1]);

		currPlot = plotHeat[2];
		currPlot[iZ] = heat[idxIn];

		currPlot = plotHeatLog[2];
		currPlot[iZ] = heat_log[idxIn];

		currPlot = plotFluence[2];
		currPlot[iZ] = fluence[idxIn];

		currPlot = plotFluenceLog[2];
		currPlot[iZ] = fluence_log[idxIn];
	}

	return;
}

// updates the slice for a given index at a given position
void mc::update_slice(const uint8_t iDim, const uint32_t idx)
{
	float* currSlice = slice[iDim]; // get pointer to slice
	float* currSliceFluence = sliceFluence[iDim];

	uint32_t idxIn, idxOut;

	if (iDim == 0) // slice me along yz
	{
		for (uint32_t iZ = 0; iZ < dims[2]; iZ++)
		{
			for (uint32_t iY = 0; iY < dims[1]; iY++)
			{
				idxIn = idx + iY * dims[0] + iZ * dims[0] * dims[1];
				idxOut = iZ + dims[2] * iY;
				currSlice[idxOut] = heat[idxIn];
				currSliceFluence[idxOut] = fluence[idxIn];
			}
		}
	}
	else if (iDim == 1) // slice me along xz
	{
		for (uint32_t iX = 0; iX < dims[0]; iX++)
		{
			for (uint32_t iZ = 0; iZ < dims[2]; iZ++)
			{
				idxIn =iX + idx * dims[0] + iZ * dims[0] * dims[1];
				idxOut = iZ + dims[2] * iX;
				currSlice[idxOut] = heat[idxIn];
				currSliceFluence[idxOut] = fluence[idxIn];
			}
		}
	}
	else if (iDim == 2) // slice me along xy
	{
		for (uint32_t iX = 0; iX < dims[0]; iX++)
		{
			for (uint32_t iY = 0; iY < dims[1]; iY++)
			{
				idxIn = iX + iY * dims[0] + idx * dims[0] * dims[1];
				idxOut = iX + dims[0] * iY;
				currSlice[idxOut] = heat[idxIn];
				currSliceFluence[idxOut] = fluence[idxIn];
			}
		}
	}
	else
	{
		printf("Invalid dimension order, must be 0, 1, or 3");
		throw "InvalidVal";
	}

	return;
}

// update logarithmic slices
void mc::update_slice_log(const uint8_t iDim, const uint32_t idx)
{
	float* currSlice = sliceLog[iDim]; // get pointer to slice
	float* currSliceFluence = sliceFluenceLog[iDim]; // get pointer to flu slice
	uint32_t idxIn, idxOut;

	if (iDim == 0) // slice me along yz
	{
		for (uint32_t iZ = 0; iZ < dims[2]; iZ++)
		{
			for (uint32_t iY = 0; iY < dims[1]; iY++)
			{
				idxOut = iZ + dims[2] * iY;
				idxIn = idx + iY * dims[0] + iZ * dims[0] * dims[1];
				currSlice[idxOut] = heat_log[idxIn];
				currSliceFluence[idxOut] = fluence_log[idxIn];
			}
		}
	}
	else if (iDim == 1) // slice me along xz
	{
		for (uint32_t iX = 0; iX < dims[0]; iX++)
		{
			for (uint32_t iZ = 0; iZ < dims[2]; iZ++)
			{
				idxOut = iZ + dims[2] * iX;
				idxIn = iX + idx * dims[0] + iZ * dims[0] * dims[1];
				currSlice[idxOut] = heat_log[idxIn];
				currSliceFluence[idxOut] = fluence_log[idxIn];
			}
		}
	}
	else if (iDim == 2) // slice me along xy
	{
		for (uint32_t iX = 0; iX < dims[0]; iX++)
		{
			for (uint32_t iY = 0; iY < dims[1]; iY++)
			{
				idxIn = iX + iY * dims[0] + idx * dims[0] * dims[1];
				idxOut = iX + dims[0] * iY;
				currSlice[idxOut] = heat_log[idxIn];
				currSliceFluence[idxOut] = fluence_log[idxIn];
			}
		}
	}
	else
	{
		printf("Invalid dimension order, must be 0, 1, or 3");
		throw "InvalidVal";
	}

	idxCross[iDim] = idx;

	return;
}

// initialize all required variables and allocate our memory
void mc::init_vars()
{	
	cudaSetDevice(sim.get_gpuID());
	cudaError_t err;
	// alloc memory for heat map and log version of it on CPU
	if (isHeatAlloc)
	{
		delete[] heat;
		delete[] heat_log;
		delete[] fluence;
		delete[] fluence_log;
		for (uint8_t iDim = 0; iDim < 3; iDim++)
		{
			delete[] slice[iDim];
			delete[] sliceLog[iDim];
			delete[] vecs[iDim];	
			delete[] plotHeat[iDim];
			delete[] plotHeatLog[iDim];
			delete[] plotFluence[iDim];
			delete[] plotFluenceLog[iDim];
		}
	}

	// allocate memory for heat and fluence on host
	heat = new float [volume.get_nElem() + 1]; // last element is container for overfly
	heat_log = new float [volume.get_nElem()];
	fluence = new float [volume.get_nElem()];
	fluence_log = new float [volume.get_nElem()];

	// grab dimensions etc from volume
	for (uint8_t iDim = 0; iDim < 3; iDim++)
	{
		dims[iDim] = volume.get_dim(iDim);
		res[iDim] = volume.get_res(iDim);
		origin[iDim] = volume.get_min(iDim);
	}

	// allocate memory for vectors and their plots
	for (uint8_t iDim = 0; iDim < 3; iDim++)
	{
		vecs[iDim] = new float [dims[iDim]];
		float* currVec = vecs[iDim];
		for (uint32_t iXYZ = 0; iXYZ < dims[iDim]; iXYZ++)
			currVec[iXYZ] = origin[iDim] + res[iDim] + ((float) iXYZ + 0.5);

		plotHeat[iDim] = new float [dims[iDim]];
		plotHeatLog[iDim] = new float [dims[iDim]];
		plotFluence[iDim] = new float [dims[iDim]];
		plotFluenceLog[iDim] = new float [dims[iDim]];
	}

	// allocate memory for slices
	slice[0] = new float [dims[1] * dims[2]]; // along x
	slice[1] = new float [dims[0] * dims[2]]; // along y
	slice[2] = new float [dims[0] * dims[1]]; // along z
	sliceLog[0] = new float [dims[1] * dims[2]]; // along x
	sliceLog[1] = new float [dims[0] * dims[2]]; // along y
	sliceLog[2] = new float [dims[0] * dims[1]]; // along z

	sliceFluence[0] = new float [dims[1] * dims[2]]; // along x
	sliceFluence[1] = new float [dims[0] * dims[2]]; // along y
	sliceFluence[2] = new float [dims[0] * dims[1]]; // along z
	sliceFluenceLog[0] = new float [dims[1] * dims[2]]; // along x
	sliceFluenceLog[1] = new float [dims[0] * dims[2]]; // along y
	sliceFluenceLog[2] = new float [dims[0] * dims[1]]; // along z

	isHeatAlloc = 1;

	for (uint64_t idx = 0; idx < volume.get_nElem(); idx++)
	{
		heat[idx] = 0;
		heat_log[idx] = 0;
		fluence[idx] = 0;
		fluence_log[idx] = 0;
	}
	heat[volume.get_nElem()] = 0; // set collection bin also to 0

	// allocate memory for our crosssection through the heat map

	// allocate matrix for heat map on GPU
	if (isHeatDevAlloc)
		cudaFree(heat_dev);

	// allocate memory on device for absorption map
	err = cudaMalloc( (void**)&heat_dev, (volume.get_nElem() + 1) * sizeof(float) );
	if (err != cudaSuccess)
	{
		printf("[mc] Could not allocate required memory for heat on card\n");
		printf("[mc] Size of requested array: nElements = %d\n", 
			volume.get_nElem());
		printf(cudaGetErrorString(err));
		printf("\n");
		throw "MemoryAllocError";
	}

	isHeatDevAlloc = 1;

	// copy zero vectors over for fluence and heat
	bool cpy1 = (cudaSuccess != cudaMemcpy(heat_dev, heat, 
		(volume.get_nElem() + 1) * sizeof(float), cudaMemcpyHostToDevice));

	if (cpy1){
		printf("Could not copy memory to card\n");	
		throw "MemoryCopyError";
	}

	return;
}

// starts the actual simulation
void mc::run_sim()
{
	cudaSetDevice(sim.get_gpuID());

	clock_t begin = clock(); // start stopwatch
	cudaError_t err; // variable used to handle CUDA errors

		// define constant simulation properties for our field	
	constArgsIn inArgs;
	for (uint8_t iDim = 0; iDim < 3; iDim++)
	{
		inArgs.dim[iDim] = volume.get_dim(iDim);
		inArgs.res[iDim] =  volume.get_res(iDim);
		inArgs.origin[iDim] = volume.get_min(iDim);
		inArgs.maxPos[iDim] = volume.get_max(iDim);
	}
	inArgs.bgMaterial = volume.get_bgMaterialId();
	inArgs.killFlag = flagKillBound;
	inArgs.nElements = inArgs.dim[0] * inArgs.dim[1] * inArgs.dim[2];

	// define optical properties of our materials
	int nMaterial = tissues.size();
	optProps* optPropHost = new optProps[nMaterial];
	for (uint8_t iMaterial = 0; iMaterial < nMaterial; iMaterial++)
	{
		// make sure that there is no zero absorption
		if (tissues[iMaterial].get_mua() < 1e-5)
			tissues[iMaterial].set_mua(1e-5);

		optPropHost[iMaterial].mu_a = tissues[iMaterial].get_mua();
		optPropHost[iMaterial].mu_s = tissues[iMaterial].get_mus();
		optPropHost[iMaterial].g = tissues[iMaterial].get_g();

		// calculated values
		optPropHost[iMaterial].albedo = tissues[iMaterial].get_albedo();
		optPropHost[iMaterial].albedo1 = 1 - optPropHost[iMaterial].albedo;
		optPropHost[iMaterial].mu_as = 1 / (optPropHost[iMaterial].mu_a + optPropHost[iMaterial].mu_s);
		optPropHost[iMaterial].g2 = optPropHost[iMaterial].g * optPropHost[iMaterial].g;
		optPropHost[iMaterial].gx2 = 2.0 * optPropHost[iMaterial].g;
	}

	// allocate memory for fiber properties
	int nFibers = fibers.size(); // number of fibers
	fiberProps* fiberPropHost = new fiberProps[nFibers];

	float totalWeight = 0;
	for (uint8_t iFiber = 0; iFiber < nFibers; iFiber++)
		totalWeight += fibers[iFiber].get_weight();
	
	float lastWeight = 0;
	for (uint8_t iFiber = 0; iFiber < nFibers; iFiber++)
	{
		// normalize fiber orientation
		const float normOr = sqrt(
			fibers[iFiber].get_orientation(0) * fibers[iFiber].get_orientation(0) +
			fibers[iFiber].get_orientation(1) * fibers[iFiber].get_orientation(1) +
			fibers[iFiber].get_orientation(2) * fibers[iFiber].get_orientation(2));

		for (uint8_t iDim = 0; iDim < 3; iDim++)
		{
			fiberPropHost[iFiber].pos[iDim] = fibers[iFiber].get_pos(iDim);
			fiberPropHost[iFiber].orientation[iDim] = 
				fibers[iFiber].get_orientation(iDim) / normOr;
		}
		fiberPropHost[iFiber].rCore = fibers[iFiber].get_rCore();
		fiberPropHost[iFiber].rCore2  = fiberPropHost[iFiber].rCore * 
			fiberPropHost[iFiber].rCore;
		fiberPropHost[iFiber].na = fibers[iFiber].get_numAp();
		
		// TODO replace this with the actual calculation depending on the medium
		fiberPropHost[iFiber].phiMax = fibers[iFiber].get_theta(1.33);
		
		// do decide which fiber emmits photon package we generate a random number
		// between 0 and 1 and project it onto a rated fiber range
		lastWeight += fibers[iFiber].get_weight() / totalWeight;
		fiberPropHost[iFiber].randMax = lastWeight;
	}
	fiberPropHost[nFibers - 1].randMax = 1.001;

	// allocate memory for tissue mapping
	uint8_t* tissueTypes_dev;
	err = cudaMalloc( (void**)&tissueTypes_dev, volume.get_nElem() * sizeof(uint8_t) );
	if (err != cudaSuccess)
	{
		printf("Could not allocate memory on card for tissue type map.\n");
		throw "cudaMallocErr";
	}

	cudaMemcpy(tissueTypes_dev, volume.get_pvolumeData(), volume.get_nElem() * sizeof(uint8_t), cudaMemcpyHostToDevice);
	if (err != cudaSuccess)
	{
		printf("Could not copy settings struct over to GPU.\n");
		throw "cudaMemcpyErr";
	}	

	// allocate memory on GPU for constant arguments and copy them over
	constArgsIn* inArgs_dev;
	err = cudaMalloc( (void**)&inArgs_dev, sizeof(constArgsIn) );
	if (err != cudaSuccess)
	{
		printf("Could not allocate memory on card for settings struct.\n");
		throw "cudaMallocErr";
	}

	cudaMemcpy(inArgs_dev, &inArgs, sizeof(constArgsIn), cudaMemcpyHostToDevice);
	if (err != cudaSuccess)
	{
		printf("Could not copy settings struct over to GPU.\n");
		throw "cudaMemcpyErr";
	}

	// allocate memory on GPU for tissue properties and copy them over
	optProps* optProp_dev;
	err = cudaMalloc( (void**)&optProp_dev, nMaterial * sizeof(optProps));
	if (err != cudaSuccess)
	{
		printf("Could not allocate memory on card for tissue struct.\n");
		throw "cudaMallocErr";
	}

	cudaMemcpy(optProp_dev, optPropHost, nMaterial * sizeof(optProps), cudaMemcpyHostToDevice);
	if (err != cudaSuccess)
	{
		printf("Could not copy tissue struct over to GPU.\n");
		throw "cudaMemcpyErr";
	}

	// alocate memory on GPU for fiber properties and copy them over
	fiberProps* fiberProp_dev;
	err = cudaMalloc( (void**)&fiberProp_dev, sizeof(fiberProps) * nFibers);
	if (err != cudaSuccess)
	{
		printf("Could not allocate memory on card for fiber struct.\n");
		throw "cudaMallocErr";
	}

	cudaMemcpy(fiberProp_dev, fiberPropHost, sizeof(fiberProps) * nFibers, cudaMemcpyHostToDevice);
	if (err != cudaSuccess)
	{
		printf("Could not copy fiber struct over to GPU.\n");
		throw "cudaMemcpyErr";
	}

 	// check occupancy before starting

  // These variables are used to convert occupancy to warps
  // cudaDeviceProp prop;
  // cudaGetDeviceProperties(&prop, sim.get_gpuID());
    
  // cudaOccupancyMaxActiveBlocksPerMultiprocessor(
  // 	&numBlocks,
  //   simPhoton,
  //   sim.get_threadsPerBlock(),
  //   0);

  int blockSize;      // The launch configurator returned block size
  int minGridSize;    // The minimum grid size needed to achieve the
          
  cudaOccupancyMaxPotentialBlockSize(
      &minGridSize,
      &blockSize,
      (void*) simPhoton,
      0,
      sim.get_nPhotonsTrue());

  int gridSize = ((sim.get_nPhotonsTrue()) + blockSize - 1) / blockSize;

	// start actual simulation
	printf("Starting actual simulation with %d blocks, each %d threads\n",
		gridSize, blockSize);

	simPhoton<<<gridSize, blockSize >>>(
		heat_dev, // output matrix into which we write our absorption
		inArgs_dev, // constant simulation parameters
		tissueTypes_dev, // defines the distribution of our tissue types
		optProp_dev, // random number 
		fiberProp_dev); // fiber properties
	cudaDeviceSynchronize();	

	err = cudaGetLastError();
	if (err != cudaSuccess){
		printf("Error occured during simulation: ");
		printf(cudaGetErrorString(err));
	}else{	

		// copy heat array back from gpu
		err = cudaMemcpy(heat, heat_dev, 
			(volume.get_nElem() + 1) * sizeof(float), cudaMemcpyDeviceToHost);
		if (err != cudaSuccess)
		{
			printf("Could not copy heat array back from GPU\n");
			throw "MemoryCopyError";
		}

		clock_t end = clock();
	
		simTime = (end - begin) / (double) CLOCKS_PER_SEC;
		//printf("Time elapsed in [s]: %f \n", simTime);
		//printf("Time per photon package in [ns]: %f \n", simTime / (double) sim.get_nPhotonsTrue() * 1e9);
		//printf("Photon packages per ms: %f \n", (double) sim.get_nPhotonsTrue() / (simTime * 1000));	
	}
	
	// free host and device memory
	cudaFree(heat_dev); // free device memorz for heat_dev
	cudaFree(fluence_dev);
	isHeatDevAlloc = 0;

	cudaFree(optProp_dev);
	cudaFree(fiberProp_dev);
	delete[] optPropHost;
	delete[] fiberPropHost;

	cudaFree(inArgs_dev);

	// for debugging: calculate and print sum of heat
	float sum = 0;
	uint32_t nanCount = 0;
	for (uint32_t iElem = 0; iElem <= volume.get_nElem(); iElem++)
	{
		if (!isnan(heat[iElem]))
			sum += heat[iElem];
		else
			nanCount++;
	}
	photonRatio = sum / ((float) sim.get_nPhotonsTrue());
	if (flagDebug)
	{
		printf("[debug] Overall photons found (0 ... 1): %.2f, nans: %d\n", 
			photonRatio, nanCount);
		printf("[debug] Percentage in collection bin (0 ... 1): %.2f\n", 
			heat[volume.get_nElem()] / sum);
	}

	// scale heat by volume and nphotons
	const float scalingFac = (float) sim.get_nPhotonsTrue() * volume.get_volume_voxel();
	uint8_t* volumeData = volume.get_pvolumeData();
	for (uint32_t iElem = 0; iElem < volume.get_nElem(); iElem++)
	{
		heat[iElem] /= scalingFac;
		const float muaTemp = tissues[volumeData[iElem]].get_mua();
		fluence[iElem] = heat[iElem] / muaTemp;
	}
	heat[volume.get_nElem()] /= scalingFac;

	calcMinMax(); // calculate minimum and maximum value in heat map
	if (flagDebug)
	{
	 	printf("[debug] maximum value in field: %f, minimum value in field: %f\n", 
	 		maxVal, minVal);
	}
	
	calcLog(); // calculate logarthmic represntation of field
	for (uint8_t iDim = 0; iDim < 3; iDim++) // generate an initial set of slices
	{
		update_slice(iDim, 0);
		update_slice_log(iDim, 0);	
	}
	
	return;
}

// claulates the minimum and the maximum value in our volume both for log and linear
// scale
void mc::calcMinMax()
{
	// get maximum and minimum value in our mip
	minVal = heat[0];
	maxVal = heat[0];
	minFluence = fluence[0];
	maxFluence = fluence[0];

	minValLog = 0; // log10(heat[0])
	maxValLog = 0; // log10(heat[0])
	minFluenceLog = 0;
	maxFluenceLog = 0;

	for (uint64_t iElement = 1; iElement < volume.get_nElem(); iElement++)
	{

		// check if we have a new maximum for our heat vector
		if (heat[iElement] > maxVal)
		{
			maxVal = heat[iElement];
			if (maxVal > 0)
			{
				maxValLog = log10(maxVal);
			}
		}

		// check if we have a new maximum for our fluence vector
		if (fluence[iElement] > maxFluence)
		{
			maxFluence = fluence[iElement];
			if (maxFluence > 0)
			{
				maxFluenceLog = log10(fluence[iElement]);
			}
		}

		// check if we have a new minimum for our heat vector
		if (heat[iElement] < minVal)
		{
			minVal = heat[iElement];
			// only update our log min if we are bigger then 0
			if (minVal > 0)
			{
				minValLog = log10(minVal);
			}
		}

		// check if we have a new minimum for our heat vector
		if (fluence[iElement] < minFluence)
		{
			minFluence = fluence[iElement];
			// only update our log min if we are bigger then 0
			if (minFluence > 0)
			{
				minFluenceLog = log10(minFluence);
			}
		}
	}

	return;
}

// converts the calculated heat map into logarthmic scale
void mc::calcLog()
{
	for (uint64_t iElement = 0; iElement < volume.get_nElem(); iElement++)
	{
		// calc logarithmic value for heat
		if (heat[iElement] > 0)
			heat_log[iElement] = log10(heat[iElement]);
		else
			heat_log[iElement] = -38; // lower end of float accuracy

		// calc logartihmic value for fluence
		if (fluence[iElement] > 0)
			fluence_log[iElement] = log10(fluence[iElement]);
		else
			fluence_log[iElement] = -38;
	}
	return;
}

void mc::run()
{
	// calc_reflectance(nWater, tissue.get_n());
	isDone = 0;
	init_vars();
	run_sim();
	isDone = 1;
	return;
}

// exports heat map as a vtk file
bool mc::exportVtk(const string filePath)
{
	bool success = 1;
	try{
		vtkwriter myWriter;

		const string title ("fluence");
		myWriter.set_title(title);
		const string type ("STRUCTURED_POINTS");
		myWriter.set_type(type);
		myWriter.set_outputPath(filePath);

		griddedData myData;

		// todo needs fixing here
		myData.data = heat;

		for (uint8_t iDim = 0; iDim < 3; iDim++)
		{
			myData.origin[iDim] = origin[iDim]; // origin in z
			myData.res[iDim] = res[iDim]; // resolution in z direction
			myData.dim[iDim] = dims[iDim];
		}

		myWriter.set_structuredPoints(&myData);
		myWriter.set_binary();

		if (flagDebug)
			printf("[debug] Writing data to file\n");

		myWriter.write();
	}
	catch(...)
	{
		printf("Could absolutely not save dataset to file\n");
	}

	return success;
}

// save fluence, heat, dimensions, resolution and origin as h5 
bool mc::exportH5(const string filePath)
{
	bool success = 1;
	try{

		H5::H5File file(filePath, H5F_ACC_TRUNC);

		// write resolutiion to file
		//printf("Write resolution to file...\n");
		const hsize_t col_dims = 3;
		H5::DataSpace mspaceRes(1, &col_dims);
		H5::DataSet resDataset = file.createDataSet(
			"dr", H5::PredType::NATIVE_FLOAT, mspaceRes);
		resDataset.write(res, H5::PredType::NATIVE_FLOAT);
		resDataset.close();

		// write origin to file
		//printf("Write origin to file...\n");
		H5::DataSpace mspaceOrigin(1, &col_dims);
		H5::DataSet originDataset = file.createDataSet(
			"origin", H5::PredType::NATIVE_FLOAT, mspaceOrigin);
		originDataset.write(origin, H5::PredType::NATIVE_FLOAT);
		originDataset.close();

		// write dimension to file
		//printf("Write dimensions to file...\n");
		H5::DataSpace mspaceDim(1, &col_dims);
		H5::DataSet dimDataset = file.createDataSet(
			"dim", H5::PredType::NATIVE_UINT, mspaceDim);
		dimDataset.write(dims, H5::PredType::NATIVE_UINT32);
		dimDataset.close();

		// write heat map to file
		const hsize_t col_data = dims[0] * dims[1] * dims[2];
		H5::DataSpace mspaceData(1, &col_data);
		H5::DataSet dataDataset = file.createDataSet(
			"heat", H5::PredType::NATIVE_FLOAT, mspaceData);
		dataDataset.write(heat, H5::PredType::NATIVE_FLOAT);
		dataDataset.close();

		// write fluence map to file
		H5::DataSpace mspaceFluence(1, &col_data);
		H5::DataSet fluenceDataset = file.createDataSet(
			"fluence", H5::PredType::NATIVE_FLOAT, mspaceFluence);
		fluenceDataset.write(fluence, H5::PredType::NATIVE_FLOAT);
		fluenceDataset.close();

		file.close();
	}
	catch(...)
	{
		printf("Could absolutely not save data to file\n");
		success = 0;
	}

	return success;
}

// save dataset as h5 to default path
bool mc::exportH5()
{
	string filePath = "/home/hofmannu/r_sync/hfoam_raw_data/Fluence/default.h5";
	bool success = exportH5(filePath);
	return success;
}

// writes all the settings of the reconstruction procedure to a file
void mc::write_settings(const string filePath)
{
	printf("Saving full simulation settings to file");
	json j;

	// simulation properties
	j["nPhotons"] = sim.get_nPhotons();
	j["killFlag"] = flagKillBound;

	// fiber properties
	// j["numAp"] = fiber->get_numAp();
	// j["dCore"] = fiber->get_dCore();
	// j["posX"] = fiber->get_pos(0);
	// j["posY"] = fiber->get_pos(1);
	// j["posZ"] = fiber->get_pos(2);
	// j["dirX"] = fiber->get_orientation(0);
	// j["dirY"] = fiber->get_orientation(1);
	// j["dirZ"] = fiber->get_orientation(2);

	std::ofstream o(filePath);
	o << j << std::endl;

	return;
}

void mc::read_settings(const string filePath)
{
	printf("Not implemented yet");
	return;
}
