/*
	Author: Urs Hofmann
	Mail: hofmannu@ethz.ch
	Date: 27.07.2021
*/


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
