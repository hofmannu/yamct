#include "optVolume.h"

// class constructor
optVolume::optVolume()
{
	// nothing here yet
}

// class destructor
optVolume::~optVolume()
{
	if (isVolumeDataAlloc)
	{
		delete[] volumeData;
		for (uint8_t iDim = 0; iDim < 3; iDim++)
		{
			delete[] vec[iDim];
			delete[] crossSection[iDim];
		}
	}
}

void optVolume::set_bgMaterialId(const uint8_t _bgMaterialId)
{
	bgMaterialId = _bgMaterialId;
	return;
}

// TODO make boundaries for checking only range from + / - radius in each dimension
void optVolume::addSphere(sphere* iSphere)
{
	float currPos[3];
	for (uint32_t ix = 0; ix < dims[0]; ix++)
	{
		currPos[0] = vec[0][ix];
		for (uint32_t iy = 0; iy < dims[1]; iy++)
		{
			currPos[1] = vec[1][iy];
			for (uint32_t iz = 0; iz < dims[2]; iz++)
			{
				currPos[2] = vec[2][iz];
				if (iSphere->isContained(&currPos[0]))
				{
					uint32_t idx = ix + iy * dims[0] + iz * dims[0] * dims[1];
					volumeData[idx] = iSphere->get_tType();
				}
			}
		}
	}

	// check if material type is higher then max annd if so, update
	if (iSphere->get_tType() > maxMaterial)
		maxMaterial = iSphere->get_tType();

	// check if material type is lower then minimum and if so, update
	if (iSphere->get_tType() < minMaterial)
		minMaterial = iSphere->get_tType();

	return;
}

// TODO make boundaries for checking only range from + / - radius in each dimension
void optVolume::addBox(box* iBox)
{
	float currPos[3];
	for (uint32_t ix = 0; ix < dims[0]; ix++)
	{
		currPos[0] = vec[0][ix];
		for (uint32_t iy = 0; iy < dims[1]; iy++)
		{
			currPos[1] = vec[1][iy];
			for (uint32_t iz = 0; iz < dims[2]; iz++)
			{
				currPos[2] = vec[2][iz];
				if (iBox->isContained(&currPos[0]))
				{
					uint32_t idx = ix + dims[0] * (iy + iz * dims[1]);
					volumeData[idx] = iBox->get_tType();
				}
			}
		}
	}

	// check if material type is higher then max and if so, update
	if (iBox->get_tType() > maxMaterial)
		maxMaterial = iBox->get_tType();

	// check if material type is higher then min value and if so, update
	if (iBox->get_tType() < minMaterial)
		minMaterial = iBox->get_tType();
	
	return;
}


// adds a tube shape to our volume
void optVolume::addTube(tube* iTube)
{
	float currPos[3];
	for (uint32_t ix = 0; ix < dims[0]; ix++)
	{
		currPos[0] = vec[0][ix];
		for (uint32_t iy = 0; iy < dims[1]; iy++)
		{
			currPos[1] = vec[1][iy];
			for (uint32_t iz = 0; iz < dims[2]; iz++)
			{
				currPos[2] = vec[2][iz];
				if (iTube->isContained(&currPos[0]))
				{
					uint32_t idx = ix + dims[0] * (iy + iz * dims[1]);
					volumeData[idx] = iTube->get_tType();
				}
			}
		}
	}

	// check if material type is higher then max and if so, update
	if (iTube->get_tType() > maxMaterial)
		maxMaterial = iTube->get_tType();

	// check if material type is higher then min value and if so, update
	if (iTube->get_tType() < minMaterial)
		minMaterial = iTube->get_tType();
	
	return;
}

// calculate size of volume along each dimension
void optVolume::calc_dims()
{
	for (uint8_t iDim = 0; iDim < 3; iDim++)
	{
		// check if upper corner value higher then lower, otherwise swap
		if (upperCorner[iDim] < lowerCorner[iDim])
		{
			printf("[optVolume] Flipping boundaries for you");
			float temp = upperCorner[iDim];
			upperCorner[iDim] = lowerCorner[iDim];
			lowerCorner[iDim] = temp;
		}

		// check if upper and lower corner have same value because then we throw an error
		if (upperCorner[iDim] == lowerCorner[iDim])
		{
			printf("[optVolume] Volume has zero size along dim %d \n", iDim);
			throw "invalidValue";
		}

		// check if resolution is valid
		if (res[iDim] == 0)
		{
			printf("[optVolume] Resolution is zero along dim %d \n", iDim);
			throw "invalidValue";
		}

		// if resolution is negative just invert 
		if (res[iDim] < 0)
		{
			res[iDim] = -res[iDim];
		}

		dims[iDim] = ceil((upperCorner[iDim] - lowerCorner[iDim]) / res[iDim]); 
	}

	nElem = dims[0] * dims[1] * dims[2];
	return;
}

// allocate volume for material index, vectors and cross sections and fill with
// bg material
void optVolume::alloc()
{
	calc_dims();

	if (isVolumeDataAlloc)
	{
		delete[] volumeData;
		delete vec[0];
		delete vec[1];
		delete vec[2];
		delete crossSection[0];
		delete crossSection[1];
		delete crossSection[2];
	}

	// generate volume data and set to background material
	volumeData = new uint8_t[nElem];
	for (uint32_t iElem = 0; iElem < nElem; iElem++)
		volumeData[iElem] = bgMaterialId;


	// allocate memory for vectors and fill them with useful information
	for (uint8_t iDim = 0; iDim < 3; iDim++)
	{
		vec[iDim] = new float [dims[iDim]];
		for (uint32_t iN = 0; iN < dims[iDim]; iN++)
		{
			vec[iDim][iN] = lowerCorner[iDim] + res[iDim] * ((float) iN);
		}
	}

	crossSection[0] = new float [dims[1] * dims[2]];
	crossSection[1] = new float [dims[0] * dims[2]];
	crossSection[2] = new float [dims[0] * dims[1]];
	
	isVolumeDataAlloc = 1;

	// reset maxMaterial value to bgMaterialId
	maxMaterial = bgMaterialId;
	minMaterial = bgMaterialId;

	return;
}

// maximum position of our volume
float optVolume::get_max(const uint8_t iDim)
{
	float maxPos = lowerCorner[iDim] + res[iDim] * ((float) dims[iDim]);
	return maxPos;
}

void optVolume::update_crossSection(const uint8_t iDim, const float crossVal)
{
	uint32_t crossIdx = (crossVal - lowerCorner[iDim]) / res[iDim];
	uint32_t idxOut, idxIn;

	if (iDim == 0) // crossIdx is x
	{
		for (uint32_t iY = 0; iY < dims[1]; iY++)
		{
			for (uint32_t iZ = 0; iZ < dims[2]; iZ++)
			{
				idxOut = iY + iZ * dims[1]; 
				idxIn = crossIdx + iY * dims[0] + iZ * dims[0] * dims[1];
				crossSection[iDim][idxOut] = volumeData[idxIn];

			}
		}
	}
	else if (iDim == 1) // crossIdx is y
	{
		for (uint32_t iX = 0; iX < dims[0]; iX++)
		{
			for (uint32_t iZ = 0; iZ < dims[2]; iZ++)
			{
				idxOut = iX + iZ * dims[0];
				idxIn = iX + crossIdx * dims[0] + iZ * dims[0] * dims[1];
				crossSection[iDim][idxOut] = volumeData[idxIn];
			}
		}
	}
	else if (iDim == 2) // crossIdx is z
	{
		for (uint32_t iX = 0; iX < dims[0]; iX++)
		{
			for (uint32_t iY = 0; iY < dims[1]; iY++)
			{
				idxOut = iX + iY * dims[0];
				idxIn = iX + iY * dims[0] + crossIdx * dims[0] * dims[1];
				crossSection[iDim][idxOut] = volumeData[idxIn];
			}
		}
	}

	return;
}


void optVolume::set_cross(const uint8_t iDim, const float crossVal)
{
	if (crossVal != cross[iDim])
	{
		cross[iDim] = crossVal;
		update_crossSection(iDim, crossVal);
	}
	return;
}

// returns the volume of each voxel in mm^3
float optVolume::get_volume_voxel()
{
	const float volume = res[0] * res[1] * res[2]; 
	return volume;
}