#include "box.h"

void box::set_tType(const uint8_t _tType)
{
	tType = _tType;
	return;
}

bool box::isContained(const float* pos)
{
	float deltaA[3];
	float deltaB[3];
	
	for (uint8_t iDim = 0; iDim < 3; iDim++)
	{
		deltaA[iDim] = pos[iDim] - cornerA[iDim];
		deltaB[iDim] = pos[iDim] - cornerB[iDim];
	}
	
	bool isContained = 0; // lets assume it is not inside

	if ((deltaA[0] <= 0) != (deltaB[0] < 0))
	{
		if ((deltaA[1] <= 0) != (deltaB[1] < 0))
		{
			if ((deltaA[2] <= 0) != (deltaB[2] < 0))
			{
				isContained = 1;
			}
		}			
	}
	
	return isContained;
}


