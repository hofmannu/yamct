// structure holding a 3d vector used for positions and velocities in x y z
__device__ struct vec3
{
	float x, y, z;
};

__device__ __inline__ float getLengthVec(const vec3 input)
{
	const float length = sqrtf(
		fmaf(input.x, input.x, fmaf(input.y, input.y, input.z * input.z)));
	return length;
}

__device__ __inline__ float getLengthVec(const vec3* input)
{
	const float length = sqrtf(
		fmaf(input->x, input->x, fmaf(input->y, input->y, input->z * input->z)));
	return length;
}

// default structure of a quaternion
__device__ struct quaternion
{
	float x, y, z, w;
};

// return conjugate of input quaternion
__device__ __inline__ static quaternion conjugate(quaternion quat)
{
  quat.x = -quat.x;
  quat.y = -quat.y;
  quat.z = -quat.z;
	// quat.w stays the same!
  return quat;
}

// multiply two quaternions
// switched to fmaf for increased float accuracy
__device__ __inline__ static quaternion mult(const quaternion& A, const quaternion& B)
{
  quaternion C;
  C.x = -A.z * B.y;
  C.x = fmaf(A.y, B.z, C.x);
  C.x = fmaf(A.x, B.w, C.x);
  C.x = fmaf(A.w, B.x, C.x);
  // C.x = A.w * B.x + A.x * B.w + A.y * B.z - A.z * B.y;

  C.y = -A.x * B.z;
  C.y = fmaf(A.z, B.x, C.y);
  C.y = fmaf(A.y, B.w, C.y);
  C.y = fmaf(A.w, B.y, C.y);
  // C.y = A.w * B.y + A.y * B.w + A.z * B.x - A.x * B.z ;

  C.z = -A.y * B.x;
  C.z = fmaf(A.z, B.w, C.z);
  C.z = fmaf(A.x, B.y, C.z);
  C.z = fmaf(A.w, B.z, C.z);
  // C.z = A.w * B.z + A.x * B.y + A.z * B.w - A.y * B.x ;

  C.w = -A.z * B.z;
  C.w = fmaf(-A.y, B.y, C.w);
  C.w = fmaf(-A.x, B.x, C.w);
  C.w = fmaf(A.w, B.w, C.w);
  // C.w = A.w * B.w - A.x * B.x - A.y * B.y - A.z * B.z;

  return C;
}

// perform a quaternion based rotation of our photon
__device__ __inline__ static void quaternion_rot(
		const vec3* rot, // vector around which we rotate (length 1)
		const float alpha, // rotation angle
		vec3* dir // direction to rotate
	)
{

	quaternion temp, quat_view, result;

	const float alpha2 = __fdividef(alpha, 2.0);
	const float sinalpha2 = __sinf(alpha2);

	// this vector is already by its definition 1 if rx, ry, rz has length 1 as well
	temp.x = rot->x * sinalpha2;
	temp.y = rot->y * sinalpha2;
	temp.z = rot->z * sinalpha2;
	temp.w = cosf(alpha2);

	// feed current position to first quaternion, already normalized as well
	quat_view.x = dir->x;
	quat_view.y = dir->y;
	quat_view.z = dir->z;
	quat_view.w = 0.0;

	result = mult(mult(temp, quat_view), conjugate(temp));
	
	dir->x = result.x; // update x output value
	dir->y = result.y; // update y output value
	dir->z = result.z; // update z output value
	
	return;
}