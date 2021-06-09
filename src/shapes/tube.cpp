#include "tube.h"

// define radius of our tube
void tube::set_iradius(const float _radius)
{
	radius[0] = _radius;
	return;
}

void tube::set_oradius(const float _radius)
{
	radius[1] = _radius;
	return;
}

// define starting position of our tube along a dimension
void tube::set_startPos(const float _startPos, const uint8_t iDim)
{
	startPos[iDim] = _startPos;
	return;
}

// define the stop position of our tube along a certain dimension
void tube::set_stopPos(const float _stopPos, const uint8_t iDim)
{
	stopPos[iDim] = _stopPos;
	return;
}

float dotProduct(const float* vec1, const float* vec2)
{
	return (vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2]);
}

float norm(const float* vec)
{
	const float length = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
	return pow(length, 0.5);
}

// def get_ref_angle(x1, x2, q): #x1 is reference pt, i.e. where angle is
//     return theta
float get_ref_angle(const float* x1, const float* x2, const float* q)
{

	float vector1 [3];
	float vector2 [3];

	for (uint8_t iDim = 0; iDim < 3; iDim++)
	{
		// vector1 = q - x1, pointing from x1 to q
		vector1[iDim] = q[iDim] - x1[iDim]; 
		// vector2 = x2 - x1, pointing from x1 to x2
		vector2[iDim] = x2[iDim] - x1[iDim];
	}

	float norm12 = norm(&vector1[0]) * norm(&vector2[0]);
	float dPro = dotProduct(&vector1[0], &vector2[0]);

//     theta = math.acos(np.dot(vector1, vector2) / (np.linalg.norm(vector1) * np.linalg.norm(vector2)))
	return acos(dPro / norm12);
}

// def get_scnd_angle(x1, x2, q):
//     vector1 = q - x2
//     vector2 = -(x2 - x1)
//     beta = math.acos(np.dot(vector1, vector2)/(np.linalg.norm(vector1)*np.linalg.norm(vector2)))
//     return beta
float get_scnd_angle(const float* x1, const float* x2, const float* q)
{

	float vector1 [3];
	float vector2 [3];

	for (uint8_t iDim = 0; iDim < 3; iDim++)
	{
		vector1[iDim] = q[iDim] - x2[iDim]; 
		vector2[iDim] = -(x2[iDim] - x1[iDim]);
	}

	float norm12 = norm(&vector1[0]) * norm(&vector2[0]);
	float dPro = dotProduct(&vector1[0], &vector2[0]);

	return acos(dPro / norm12);
}

// returns the orthogonal distance between a line (defined by x1 / x2) and a point
// # returns the closest distance between an infinite line (x1 to x2) and a point
// def get_orthog_distance(x1, x2, q):
//     r = np.linalg.norm(q - x1) * math.sin(get_ref_angle(x1, x2, q))
//     return r
float get_orthog_dist(const float* x1, const float* x2, const float* q)
{
	float distVec[3];
	for (uint8_t iDim = 0; iDim < 3; iDim++)
	{
		distVec[iDim] = q[iDim] - x1[iDim]; 
	}
	float dist = norm(&distVec[0]);

	float angle = get_ref_angle(x1, x2, q);

	return fabs(dist * sin(angle));
}


// returns if current position is contained within tube
// def is_in_tube(q, x1, x2, ri, ro):
//     r_orthog = get_orthog_distance(x1, x2, q)
//     angle1 = get_ref_angle(x1, x2, q)
//     angle2 = get_scnd_angle(x1, x2, q)
    
//     if (angle1 <= np.pi/2) and (angle2 <= np.pi/2):
// #         if r_orthog < ri:
// #             return 0
//         if r_orthog >= ri and r_orthog <= ro:
//             # print("transcending barrier")    
//             return 1
//         else:
//             return 0
//     else:
//         return 0
// TODO make more efficient
bool tube::isContained(const float* pos)
{
	bool isContained = 0;
	float r_orthog = get_orthog_dist(&startPos[0], &stopPos[0], pos);
	float angle1 = get_ref_angle(&startPos[0], &stopPos[0], pos);
	float angle2 = get_scnd_angle(&startPos[0], &stopPos[0], pos);

	if ((angle1 <= (M_PI / 2)) && (angle2 <= (M_PI / 2)))
	{
		if (radius[0] > 0)
		{
			if ((r_orthog >= radius[0]) && (r_orthog <= radius[1]))
			{
				isContained = 1;
			}
		}
		else
		{
			if (r_orthog <= radius[1])
			{
				isContained = 1;
			}
		}
	}

	return isContained;
}

