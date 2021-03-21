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

// generate a new sphere and add it to our already existing list
void optVolume::new_sphere()
{
	sphere newSphere;
	spheres.push_back(newSphere);
	return;
}

void optVolume::delete_sphere(const uint8_t iSphere)
{
	spheres.erase(spheres.begin() + iSphere);
	return;
}

// TODO make boundaries for checking only range from + / - radius in each dimension
void optVolume::add_sphere(sphere* iSphere)
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

void optVolume::set_sphere(const uint8_t iSphere, const sphere sphereObj)
{
	spheres[iSphere] = sphereObj;
	return;
}

void optVolume::new_box()
{
	box newBox;
	boxes.push_back(newBox);
	return;
}

void optVolume::delete_box(const uint8_t iBox)
{
	boxes.erase(boxes.begin() + iBox);
	return;
}

// TODO make boundaries for checking only range from + / - radius in each dimension
void optVolume::add_box(box* iBox)
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

void optVolume::set_box(const uint8_t iBox, const box boxObj)
{
	boxes[iBox] = boxObj;
	return;
}

void optVolume::new_tube()
{
	tube newTube;
	tubes.push_back(newTube);
	return;
}

void optVolume::delete_tube(const uint8_t iTube)
{
	tubes.erase(tubes.begin() + iTube);
	return;
}

// adds a tube shape to our volume
void optVolume::add_tube(tube* iTube)
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

void optVolume::set_tube(const uint8_t iTube, const tube tubeObj)
{
	tubes[iTube] = tubeObj;
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

void optVolume::generate_volume()
{
	for (uint8_t iPriority = 0; iPriority < 255; iPriority++)
	{
		// go through all shapes and check if they match this priority, if so push them to
		// volume class
		for (uint8_t iSphere = 0; iSphere < spheres.size(); iSphere++)
		{
			if (spheres[iSphere].get_priority() == iPriority)
			{
				printf("Adding a sphere\n");
				add_sphere(&spheres[iSphere]);
			}
		}

		for (uint8_t iBox = 0; iBox < boxes.size(); iBox++)
		{
			if (boxes[iBox].get_priority() == iPriority)
			{
				printf("Adding a box\n");
				add_box(&boxes[iBox]);
			}
		}

		for (uint8_t iTube = 0; iTube < tubes.size(); iTube++)
		{
			if (tubes[iTube].get_priority() == iPriority)
			{
				printf("Adding a tube\n");
				add_tube(&tubes[iTube]);
			}
		}

	}
	return;
}

// return all geometrioes as a json file
json optVolume::get_json()
{
	json jv;

	// write all important field properties to file
	jv["fieldProps"] = {
		{"resX", res[0]},
		{"resY", res[1]},
		{"resZ", res[2]},
		{"minX", lowerCorner[0]},
		{"minY", lowerCorner[1]},
		{"minZ", lowerCorner[2]},
		{"maxX", upperCorner[0]},
		{"maxY", upperCorner[1]},
		{"maxZ", upperCorner[2]},
		{"bgMaterialId", bgMaterialId}
	};

	jv["nSpheres"] = spheres.size();
	jv["nTubes"] = tubes.size();
	jv["nBoxes"] = boxes.size();

	// add all spheres to out config
	for (uint8_t iSphere = 0; iSphere < spheres.size(); iSphere++)
	{
		char labelSphere[80];
		sprintf(labelSphere, "sphere%d", iSphere);
		jv[labelSphere] = {
			{"r", spheres[iSphere].get_radius()},
			{"centerX", spheres[iSphere].get_center(0)},
			{"centerY", spheres[iSphere].get_center(1)},
			{"centerZ", spheres[iSphere].get_center(2)},
			{"priority", spheres[iSphere].get_priority()},
			{"material", spheres[iSphere].get_tType()}
		};
	}

	for (uint8_t iTube = 0; iTube < tubes.size(); iTube++)
	{
		char labelTube[80];
		sprintf(labelTube, "tube%d", iTube);
		jv[labelTube] = {
			{"r1", tubes[iTube].get_iradius()},
			{"r2", tubes[iTube].get_oradius()},
			{"startPosX", tubes[iTube].get_startPos(0)},
			{"startPosY", tubes[iTube].get_startPos(1)},
			{"startPosZ", tubes[iTube].get_startPos(2)},
			{"stopPosX", tubes[iTube].get_stopPos(0)},
			{"stopPosY", tubes[iTube].get_stopPos(1)},
			{"stopPosZ", tubes[iTube].get_stopPos(2)},
			{"priority",tubes[iTube].get_priority() },
			{"material",tubes[iTube].get_tType() }
		};
	}

	for (uint8_t iBox = 0; iBox < boxes.size(); iBox++)
	{
		char labelBox[80];
		sprintf(labelBox, "box%d", iBox);
		jv[labelBox] = {
			{"cornerAX", boxes[iBox].get_cornerA(0)},
			{"cornerAY", boxes[iBox].get_cornerA(1)},
			{"cornerAZ", boxes[iBox].get_cornerA(2)},
			{"cornerBX", boxes[iBox].get_cornerB(0)},
			{"cornerBY", boxes[iBox].get_cornerB(1)},
			{"cornerBZ", boxes[iBox].get_cornerB(2)},
			{"priority", boxes[iBox].get_priority()},
			{"material", boxes[iBox].get_tType()}		
		};
	}

	return jv;
}

// import all settings for the geometry from out settings file
void optVolume::read_json(const json jOpt)
{

	uint8_t nSpheres = jOpt["nSpheres"];
	uint8_t nTubes = jOpt["nTubes"];
	uint8_t nBoxes = jOpt["nBoxes"];


	res[0] = jOpt["fieldProps"]["resX"];
	res[1] = jOpt["fieldProps"]["resY"];
	res[2] = jOpt["fieldProps"]["resZ"];
	lowerCorner[0] = jOpt["fieldProps"]["minX"];
	lowerCorner[1] = jOpt["fieldProps"]["minY"]; 
	lowerCorner[2] = jOpt["fieldProps"]["minZ"];
	upperCorner[0] = jOpt["fieldProps"]["maxX"];
	upperCorner[1] = jOpt["fieldProps"]["maxY"];
	upperCorner[2] = jOpt["fieldProps"]["maxZ"];
	bgMaterialId = jOpt["fieldProps"]["bgMaterialId"];
	
	// add all spheres to out config
	spheres.clear();
	for (uint8_t iSphere = 0; iSphere < nSpheres; iSphere++)
	{
		sphere currSphere;
		char labelSphere[80];
		sprintf(labelSphere, "sphere%d", iSphere);
		currSphere.set_radius(jOpt[labelSphere]["r"]);
		currSphere.set_center(jOpt[labelSphere]["centerX"], 0);
		currSphere.set_center(jOpt[labelSphere]["centerY"], 1);
		currSphere.set_center(jOpt[labelSphere]["centerZ"], 2);
		currSphere.set_priority(jOpt[labelSphere]["priority"]);
		currSphere.set_tType(jOpt[labelSphere]["material"]);
		spheres.push_back(currSphere);
	}

	tubes.clear();
	for (uint8_t iTube = 0; iTube < nTubes; iTube++)
	{
		tube currTube;
		char labelTube[80];
		sprintf(labelTube, "tube%d", iTube);

		currTube.set_iradius(jOpt[labelTube]["r1"]);
		currTube.set_oradius(jOpt[labelTube]["r2"]);
		currTube.set_startPos(jOpt[labelTube]["startPosX"], 0);
		currTube.set_startPos(jOpt[labelTube]["startPosY"], 1);
		currTube.set_startPos(jOpt[labelTube]["startPosZ"], 2);
		currTube.set_stopPos(jOpt[labelTube]["stopPosX"], 0);
		currTube.set_stopPos(jOpt[labelTube]["stopPosY"], 1);
		currTube.set_stopPos(jOpt[labelTube]["stopPosZ"], 2);
		currTube.set_priority(jOpt[labelTube]["priority"]);
		currTube.set_tType(jOpt[labelTube]["material"]);
		tubes.push_back(currTube);
	}

	boxes.clear();
	for (uint8_t iBox = 0; iBox < nBoxes; iBox++)
	{
		box currBox;
		char labelBox[80];
		sprintf(labelBox, "box%d", iBox);
		
		currBox.set_cornerA(jOpt[labelBox]["cornerAX"], 0);
		currBox.set_cornerA(jOpt[labelBox]["cornerAY"], 1);
		currBox.set_cornerA(jOpt[labelBox]["cornerAZ"], 2);
		currBox.set_cornerB(jOpt[labelBox]["cornerBX"], 0);
		currBox.set_cornerB(jOpt[labelBox]["cornerBY"], 1);
		currBox.set_cornerB(jOpt[labelBox]["cornerBZ"], 2);
		currBox.set_priority(jOpt[labelBox]["priority"]);
		currBox.set_tType(jOpt[labelBox]["material"]);

		boxes.push_back(currBox);
	}

	return;
}