#include "interface.cuh"

// class constructor
interface::interface()
{
	simprop = sim.get_psim(); // get pointer to simulation properties
	volume = sim.get_pvolume(); // get pointer to volumetric representation
	tissueTypes = sim.get_ptissues(); // get pointer to tissue types
	fibers = sim.get_pfibers();

	// update information about GPUs
	ScanGpus();

}

// class destructor
interface::~interface()
{

}

// check available GPUs 
void interface::ScanGpus()
{
	// empty vector from old stuff

	// while constructing our interface check how many GPUs are present
	cudaGetDeviceCount(&deviceCount);
	deviceNames.clear();
	for (uint8_t iGPU = 0; iGPU < deviceCount; iGPU++)
	{
		cudaDeviceProp deviceProperties;
	  cudaGetDeviceProperties(&deviceProperties, iGPU);
	 	deviceNames.push_back(deviceProperties);
	}
	return;
}

// displays a small help marker next to the text
static void HelpMarker(const char* desc)
{
	ImGui::TextDisabled("(?)");
	if (ImGui::IsItemHovered())
	{
		ImGui::BeginTooltip();
		ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
		ImGui::TextUnformatted(desc);
		ImGui::PopTextWrapPos();
		ImGui::EndTooltip();
	}
	return;
}

void interface::InitWindow(int *argcp, char**argv)
{

	if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER | SDL_INIT_GAMECONTROLLER) != 0)
	{
	  printf("Error: %s\n", SDL_GetError());
		return;
	}
	// main_display_function goes somewhere here
	const char* glsl_version = "#version 140";
	// to find out which glsl version you are using, run glxinfo from terminal
	// and look for "OpenGL shading language version string"
	// https://en.wikipedia.org/wiki/OpenGL_Shading_Language

	SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, 0);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 0);
	
	SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
	SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
	SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);

	SDL_WindowFlags window_flags = (SDL_WindowFlags)(SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI);
	SDL_Window* window = SDL_CreateWindow(windowTitle, 
		SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 1900, 1080, window_flags);
	SDL_GLContext gl_context = SDL_GL_CreateContext(window);
	SDL_GL_MakeCurrent(window, gl_context);
	SDL_GL_SetSwapInterval(1); // Enable vsync

	bool err = glewInit() != GLEW_OK;
	if (err)
	{
	  throw "Failed to initialize OpenGL loader!";
		return;
	}

	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO(); (void)io;

	ImGui::StyleColorsDark();
	ImGui_ImplSDL2_InitForOpenGL(window, gl_context);

	ImGui_ImplOpenGL3_Init(glsl_version);
	bool done = false;
	while (!done)
	{
		SDL_Event event;
		while (SDL_PollEvent(&event))
		{
			ImGui_ImplSDL2_ProcessEvent(&event);
			if (event.type == SDL_QUIT)
				done = true;
		}
		// Start the Dear ImGui frame
		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplSDL2_NewFrame(window);
		ImGui::NewFrame();
		MainDisplayCode();
		// Rendering
		ImGui::Render();
		
		glViewport(0, 0, (int)io.DisplaySize.x, (int)io.DisplaySize.y);
		glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
		glClear(GL_COLOR_BUFFER_BIT);
		//glUseProgram(0); // You may want this if using this code in an OpenGL 3+ context where shaders may be bound
		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
		SDL_GL_SwapWindow(window);
	}

	// Cleanup
	ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplSDL2_Shutdown();
	ImGui::DestroyContext();

	SDL_GL_DeleteContext(gl_context);
 	SDL_DestroyWindow(window);
 	SDL_Quit();
	return;
}


int getSPcores(cudaDeviceProp devProp)
{  
    int cores = 0;
    int mp = devProp.multiProcessorCount;
    switch (devProp.major){
     case 2: // Fermi
      if (devProp.minor == 1) cores = mp * 48;
      else cores = mp * 32;
      break;
     case 3: // Kepler
      cores = mp * 192;
      break;
     case 5: // Maxwell
      cores = mp * 128;
      break;
     case 6: // Pascal
      if ((devProp.minor == 1) || (devProp.minor == 2)) cores = mp * 128;
      else if (devProp.minor == 0) cores = mp * 64;
      else printf("Unknown device type\n");
      break;
     case 7: // Volta and Turing
      if ((devProp.minor == 0) || (devProp.minor == 5)) cores = mp * 64;
      else printf("Unknown device type\n");
      break;
     case 8: // Ampere
      if (devProp.minor == 0) cores = mp * 64;
      else if (devProp.minor == 6) cores = mp * 128;
      else printf("Unknown device type\n");
      break;
     default:
      printf("Unknown device type\n"); 
      break;
      }
    return cores;
}

// interfacing for important properties like fiber stuff, field size etc
void interface::Properties()
{
	ImGui::Begin("Simulation", &show_properties_window);
	
	ImGui::Columns(1);
	
	ImGui::Separator();
	ImGui::Text("Simulation properties");
	int nPhotons = simprop->get_nPhotons();
	ImGui::InputInt("Number of photons", &nPhotons); ImGui::SameLine();
	HelpMarker("Number of photon packages simulated using Monte Carlo.s");
	simprop->set_nPhotons(nPhotons);
	ImGui::Text("True number of simulated photons: %d", simprop->get_nPhotonsTrue());
	ImGui::Text("Number of blocks: %d, threads per block: %d", 
		simprop->get_nBlocks(), simprop->get_threadsPerBlock());

	ImGui::Checkbox("Kill at boundary", sim.get_pflagKillBound());

	ImGui::Separator();

	ImGui::Columns(4);
	ImGui::Text("ID"); ImGui::NextColumn();
	ImGui::Text("Name"); ImGui::NextColumn();
	ImGui::Text("Cores"); ImGui::NextColumn();
	ImGui::Text("Ram"); ImGui::NextColumn();
	for (uint8_t iGPU = 0; iGPU < deviceCount; iGPU++)
	{
		cudaDeviceProp currProp = deviceNames[iGPU];
		ImGui::Text("%d", iGPU); ImGui::NextColumn();
		ImGui::Text("%s", currProp.name); ImGui::NextColumn();
		ImGui::Text("%d", getSPcores(currProp)); ImGui::NextColumn();
		float memGb = ((float) currProp.totalGlobalMem) / 1024.0 / 1024.0 / 1024.0;
		ImGui::Text("%.1f", memGb); ImGui::NextColumn();
	}
	int gpuID = simprop->get_gpuID();

	ImGui::Columns(2);
	ImGui::InputInt("GPU to use", &gpuID);
	ImGui::NextColumn();
	if (ImGui::Button("Rescan"))
	{
		ScanGpus();
	}
	ImGui::NextColumn();

	simprop->set_gpuID(gpuID);


	// check if everything is ready for our simulation
	uint8_t maxMaterial = volume->get_maxMaterial();
	uint8_t nMaterials = tissueTypes->size();
	ImGui::Separator();
	ImGui::Columns(2);

	// is a valid gpu present
	ImGui::Text("Is CUDA capable GPU present and selected?");
	ImGui::NextColumn();
	bool isGpuOk = (deviceCount > 0) && (simprop->get_gpuID() < deviceCount);
	if (isGpuOk)
		ImGui::Text("y");
	else
		ImGui::Text("n");
	ImGui::NextColumn();

	// is our volume generated 
	ImGui::Text("Is volume generated?");
	ImGui::NextColumn();
	if (is_volume_generated)
		ImGui::Text("y");
	else
		ImGui::Text("n");
	ImGui::NextColumn();

	// is every required material defined
	bool is_materials_defined = (maxMaterial < nMaterials);
	ImGui::Text("Required materials defined?");
	ImGui::NextColumn();
	if (is_materials_defined)
		ImGui::Text("y");
	else
		ImGui::Text("n");
	ImGui::NextColumn();

	ImGui::Columns(1);

	bool is_all_valid = (is_volume_generated && is_materials_defined &&
		isGpuOk);

	// if something not ready yet, disable button
	if (!is_all_valid)
	{
		ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
    ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
	}

	if (ImGui::Button("Run"))
	{
		sim.run();
		is_output_defined = 1;
			
		// ImGui::OpenPopup("Ongoing simulation");
		// if (ImGui::BeginPopupModal("Ongoing simulation"))
		// {
		// 	ImGui::Text("Progress bar missing");

		// 	// run simulation

		// 	std::thread tsim (&mc::run, &sim);


		// 	if (sim.get_isDone())
		// 	{
		// 		tsim.join();
	 //    	ImGui::CloseCurrentPopup();
		// 	}

		//   ImGui::EndPopup();
		// }
	}

	if (!is_all_valid)
	{
		ImGui::PopItemFlag();
    	ImGui::PopStyleVar();
	}

	ImGui::Separator();
	// option to save and load settings to a file

	static char filePath [64];
	ImGui::InputText("Settings path", filePath, 64);
	string filePathString = filePath;

	ImGui::Columns(2);

	if (ImGui::Button("Save settings"))
		sim.write_settings(filePath);
	
	ImGui::NextColumn();
	if (ImGui::Button("Load settings"))
	{
		// if path is pointing to a file, do the reading, otherwise throw a warning 
		if (boost::filesystem::exists(filePath))
		{
			sim.read_settings(filePath);
		}
		else
		{
			printf("path is not pointing to a file, gonna ignore request\n");
			ImGui::OpenPopup("Path not pointing to file");	
		}
	}

	if (ImGui::BeginPopupModal("Path not pointing to file"))
	{
		ImGui::Text("Path is not pointing to a file");
	    if (ImGui::Button("Close"))
	    	ImGui::CloseCurrentPopup();
	    ImGui::EndPopup();
	}

	ImGui::End();
	return;
}

// allows the definition of our illumination properties
void interface::Illumination()
{

	ImGui::Begin("Illumination", &show_illumination);
	
	if (ImGui::Button("Add fiber"))
	{
		fiberProperties newFiber;
		fibers->push_back(newFiber);
	}

	// amek editable fields for all the fiber properties
	for (uint8_t iFiber = 0; iFiber < fibers->size(); iFiber++)
	{
		ImGui::PushID(iFiber);
		ImGui::Columns(2);
		ImGui::Text("Fiber ID: %d", iFiber); ImGui::NextColumn();
		if (ImGui::Button("x"))
		{
			fibers->erase(fibers->begin() + iFiber);
		}
		else
		{
			ImGui::Columns(1);
			// fiber properties
			fiberProperties currFiber = fibers->at(iFiber);
			ImGui::InputFloat("Numerical aperture", currFiber.get_pnumAp());
			ImGui::InputFloat("Core diameter [mm]", currFiber.get_pdCore());
			ImGui::InputFloat3("Position x/y/z [mm]", currFiber.get_ppos());
			ImGui::InputFloat3("Orientation x/y/z [mm]", currFiber.get_porientation());	
			ImGui::InputFloat("Weight [1]", currFiber.get_pweight());
			fibers->at(iFiber) = currFiber;
		}
		ImGui::PopID();
		ImGui::Separator();
	}
	ImGui::End();
	return;
}

// c implementation of matlabs imagesc
void interface::ImImagesc(
	const double* data, const uint64_t sizex, const uint64_t sizey, 
	GLuint* out_texture, const color_mapper myCMap)
{
	
	glDeleteTextures(1, out_texture);

	// Create an OpenGL texture identifier
	GLuint image_texture;
	glGenTextures(1, &image_texture);
	glBindTexture(GL_TEXTURE_2D, image_texture);

	// setup filtering parameters for display
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
			
	// use color transfer function to convert from float to rgba
	unsigned char* data_conv = new unsigned char[4 * sizex * sizey];
	myCMap.convert_to_map(data, sizex * sizey, data_conv);
	glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, sizex, sizey, 0, GL_RGBA, GL_UNSIGNED_BYTE, 
		data_conv);

	// give pointer back to main program
	*out_texture = image_texture;
	delete[] data_conv; // free memory for temporary array
	return;
}


void interface::ImImagesc(
	const float* data, const uint64_t sizex, const uint64_t sizey, 
	GLuint* out_texture, const color_mapper myCMap)
{
	
	glDeleteTextures(1, out_texture);

	// Create an OpenGL texture identifier
	GLuint image_texture;
	glGenTextures(1, &image_texture);
	glBindTexture(GL_TEXTURE_2D, image_texture);

	// setup filtering parameters for display
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
			
	// use color transfer function to convert from float to rgba
	unsigned char* data_conv = new unsigned char[4 * sizex * sizey];
	myCMap.convert_to_map(data, sizex * sizey, data_conv);
	glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, sizex, sizey, 0, GL_RGBA, GL_UNSIGNED_BYTE, 
		data_conv);

	// give pointer back to main program
	*out_texture = image_texture;
	delete[] data_conv; // free memory for temporary array
	return;
}
// allows user to define tissue types with different optical properties
void interface::TissueProperties()
{
	ImGui::Begin("Optical properties", &show_tissue_properties);
	if (ImGui::Button("Add Type"))
	{
		optProperties newTissue;
		tissueTypes->push_back(newTissue);
	}

	ImGui::Columns(2);
	for (uint8_t iTissue = 0; iTissue < tissueTypes->size(); iTissue++)
	{

		ImGui::PushID(iTissue);
		ImGui::Text("Tissue ID: %d", iTissue); ImGui::NextColumn();
		if (ImGui::Button("x"))
		{
			tissueTypes->erase(tissueTypes->begin() + iTissue);
		}
		else
		{
			optProperties currTissue = tissueTypes->at(iTissue);
			ImGui::NextColumn();
			ImGui::InputFloat("Abs coeff [1/mm]", currTissue.get_pmua()); 
			ImGui::NextColumn();
			ImGui::InputFloat("Scatt coeff [1/mm]", currTissue.get_pmus()); 
			ImGui::NextColumn();
			ImGui::InputFloat("Anisotropy", currTissue.get_pg()); 
			ImGui::NextColumn();
			ImGui::InputFloat("Refractive index", currTissue.get_pn());
			ImGui::NextColumn();
			tissueTypes->at(iTissue) = currTissue;
		}
		ImGui::PopID();
		ImGui::Separator();
	}
	ImGui::Columns(1);
	ImGui::End();
	return;
}

// field generator: allows user to define different geometries of absorbing or scattering
// material on top of our bacground material
// TODO include prioritizing into options
void interface::FieldGenerator()
{
	ImGui::Begin("Field generator", &show_field_generator);

	ImGui::Text("Volume definition");
	// x range, y range, z range of simulation
	ImGui::InputFloat3("Resolution x/y/z [mm]", volume->get_pres());
	ImGui::InputFloat3("Minimum value x/y/z [mm]", volume->get_plowerCorner());
	ImGui::InputFloat3("Maximum value x/y/z [mm]", volume->get_pupperCorner());
	int bgMaterialId = volume->get_bgMaterialId();
	ImGui::InputInt("Background material", &bgMaterialId);
	volume->set_bgMaterialId(bgMaterialId);	

	// background material

	ImGui::Columns(2);
	ImGui::Separator();
	ImGui::Text("Sphere");
	ImGui::NextColumn();
	if (ImGui::Button("Add Sphere"))
	{
		volume->new_sphere();
	}
	ImGui::NextColumn();

	for (uint32_t iSphere = 0; iSphere < volume->get_nsphere(); iSphere++)
	{
		sphere* currSphere = volume->get_psphere(iSphere);
		ImGui::PushID(iSphere);
		ImGui::InputFloat("Radius [mm]", currSphere->get_pradius());
		ImGui::NextColumn();
		int tType = currSphere->get_tType();
		ImGui::InputInt("Tissue type", &tType);
		currSphere->set_tType(tType);
		ImGui::NextColumn();
		ImGui::InputFloat3("Center x/y/z [mm]", currSphere->get_pcenter());
		ImGui::NextColumn();
		if (ImGui::Button("x"))
		{
			volume->delete_sphere(iSphere);
		}
		ImGui::PopID();
		ImGui::NextColumn();
	}

	ImGui::Separator();
	ImGui::Text("Box");
	ImGui::NextColumn();
	if (ImGui::Button("Add Box"))
	{
		volume->new_box();
	}
	ImGui::NextColumn();

	for (uint32_t iBox = 0; iBox < volume->get_nbox(); iBox++)
	{
		box* currBox = volume->get_pbox(iBox);
		ImGui::PushID(iBox + 255);
		ImGui::InputFloat3("Corner A x/y/z [mm]", currBox->get_pcornerA());
		ImGui::NextColumn();
		ImGui::InputFloat3("Corner B x/y/z [mm]", currBox->get_pcornerB());
		ImGui::NextColumn();
		int tType = currBox->get_tType();
		ImGui::InputInt("Tissue type", &tType);
		currBox->set_tType(tType);
		ImGui::NextColumn();
		if (ImGui::Button("x"))
		{
			volume->delete_box(iBox);
		}
		ImGui::PopID();
		ImGui::NextColumn();
	}

	ImGui::Separator();
	ImGui::Text("Cylinder");
	ImGui::NextColumn();
	if (ImGui::Button("Add Cylinder"))
	{
		
	}
	ImGui::NextColumn();

	ImGui::Separator();
	ImGui::Text("Tubes");
	ImGui::NextColumn();
	if (ImGui::Button("Add Tube"))
	{
		volume->new_tube();
	}
	ImGui::NextColumn();

	for (uint32_t iTube = 0; iTube < volume->get_ntube(); iTube++)
	{
		tube* currTube = volume->get_ptube(iTube);
		ImGui::Columns(2);
		ImGui::PushID(iTube + 255 * 2);
		ImGui::InputFloat3("Start pos A x/y/z [mm]", currTube->get_pstartPos());
		ImGui::NextColumn();
		ImGui::InputFloat3("Stop pos B x/y/z [mm]", currTube->get_pstopPos());
		ImGui::NextColumn();
		ImGui::InputFloat2("Radius (inner, outer) [mm]", currTube->get_pradius());
		ImGui::NextColumn();
		int tType = currTube->get_tType();
		ImGui::InputInt("Tissue type", &tType);
		currTube->set_tType(tType);
		ImGui::SameLine();
		if (ImGui::Button("x"))
		{
			volume->delete_tube(iTube);
		}
		ImGui::PopID();
		ImGui::NextColumn();
	}

	ImGui::Columns(1);

	// generate volumetric representation
	if (ImGui::Button("Generate volume"))
	{
		// allocate memory for volume here
		volume->alloc();
		volume->generate_volume();


		// update crosssections through our volume
		for (uint8_t iDim = 0; iDim < 3; iDim++)
			volume->update_crossSection(iDim, volume->get_cross(iDim));

		// once done, set flag to high
		is_volume_generated = 1;

		// push max and min id as limits over to preview
		material_mapper.set_maxVal(volume->get_maxMaterial());
		material_mapper.set_minVal(volume->get_minMaterial());
	}

	// if volume was generated make a small preview here
	if (is_volume_generated)
	{
		ImGui::Text("Dimension of generated volume: %d x %d x %d", 
			volume->get_dim(0), volume->get_dim(1), volume->get_dim(2));
		ImGui::Text("Min material: %d, max material %d", 
			volume->get_minMaterial(), volume->get_maxMaterial());
		
		int width = 550;
		int height = 550;

		ImImagesc(volume->get_pcrossSection(0),
	 		volume->get_dim(1), volume->get_dim(2), &material_texture, material_mapper);
		ImGui::Image((void*)(intptr_t)material_texture, ImVec2(width, height));

		// cross section along depth

		// cross section along width

		// slider for depth goes here
		float crossX = volume->get_cross(0);
		ImGui::SliderFloat("x crossection", &crossX, 
			volume->get_min(0), volume->get_max(0), "%.2f");
		volume->set_cross(0, crossX);

		// ImGui::SliderFloat("MinVal", material_mapper.get_pminVal(), 0, 255, "%.1f");
		// ImGui::SliderFloat("MaxVal", material_mapper.get_pmaxVal(), 0, 255, "%.1f");
		ImGui::ColorEdit4("Min color", material_mapper.get_pminCol(), ImGuiColorEditFlags_Float);
		ImGui::ColorEdit4("Max color", material_mapper.get_pmaxCol(), ImGuiColorEditFlags_Float);

		// float crossY = volume.get_cross(1);
		// ImGui::SliderFloat("y crosssection", volume.get_pcross(1), 
		// 	volume.get_min(1), volume.get_max(1), "%.2f");
		// volume.set_cross(1, crossY);

		// slider for width goes here
	}

	ImGui::End();

	return;
}

void interface::Result()
{
	ImGui::Begin("Result", &show_results);

	// make imagesc here with plot below
	int overallHeight = 800;
	float rangeX = volume->get_range(0);
	float rangeY = volume->get_range(1);
	float rangeZ = volume->get_range(2);

	float heightVal = rangeX + rangeZ;
	float widthVal = rangeY;

	int width = (float) overallHeight * widthVal / heightVal;
	int height1 = (float) overallHeight / heightVal * rangeZ; // (float) width / field->get_rMax() * field->get_zExtend();
	int height2 = (float) overallHeight / heightVal * rangeX;
	// ImGui::Checkbox("Logarthmic scale", &flagLogPlot);

	ImGui::Columns(2);
	ImGui::Text("Simulation time: %.1f s", sim.get_simTime());
	ImGui::NextColumn();
	ImGui::Text("Photon ratio: %.1f percent", sim.get_photonRatio() * 100);
	ImGui::NextColumn();
	ImGui::Checkbox("Log scale?", &flagLog);
	ImGui::NextColumn();
	ImGui::Checkbox("Fluence?", &flagFluence);
	ImGui::NextColumn();
	ImGui::Columns(1);

	if (is_output_defined) // use logarithmic plot in this case
	{
		ImGui::Separator();

		float* sliceX = sim.get_slice(0, xPos, flagLog, flagFluence);
		float* sliceY = sim.get_slice(1, yPos, flagLog, flagFluence);

		ImImagesc(sliceX,
			sim.get_dim(1), sim.get_dim(2), &fluence_texture_x, fluence_mapper);
		ImGui::Image((void*)(intptr_t)fluence_texture_x, ImVec2(width, height1));
		ImImagesc(sliceY,
			sim.get_dim(2), sim.get_dim(0), &fluence_texture_y, fluence_mapper);
		ImGui::Image((void*)(intptr_t)fluence_texture_y, ImVec2(width, height2));

		float maxVal, minVal;
		if (flagFluence)
		{
			maxVal = (flagLog == 1) ? sim.get_maxFluenceLog() : sim.get_maxFluence();
			minVal = (flagLog == 1) ? sim.get_minFluenceLog() : sim.get_minFluence();
		}
		else
		{
			maxVal = (flagLog == 1) ? sim.get_maxValLog() : sim.get_maxVal();
			minVal = (flagLog == 1) ? sim.get_minValLog() : sim.get_minVal();
		}

		ImGui::SliderFloat("MinVal", fluence_mapper.get_pminVal(), 
			minVal, maxVal, "%.9f");
		ImGui::SliderFloat("MaxVal", fluence_mapper.get_pmaxVal(), 
			minVal, maxVal, "%.9f");

		ImGui::ColorEdit4("Min color", fluence_mapper.get_pminCol(), ImGuiColorEditFlags_Float);
		ImGui::ColorEdit4("Max color", fluence_mapper.get_pmaxCol(), ImGuiColorEditFlags_Float);
		
		ImGui::SliderFloat("x pos", &xPos, sim.get_minPos(0), sim.get_maxPos(0), "%.2f");
		ImGui::SliderFloat("y pos", &yPos, sim.get_minPos(1), sim.get_maxPos(1), "%.2f");
		ImGui::SliderFloat("z pos", &zPos, sim.get_minPos(2), sim.get_maxPos(2), "%.2f");

		ImGui::Separator();	
		float posXYZ [3] = {xPos, yPos, zPos};

		// // plot axial crossection
		ImGui::PlotConfig conf;
		float* plotVec = sim.get_plot(0, &posXYZ[0], flagLog, flagFluence);
		float* plotVecF = new float [sim.get_dim(0)];
		for (uint32_t iElement = 0; iElement < sim.get_dim(0); iElement++)
			plotVecF[iElement] = (float) plotVec[iElement];

		conf.values.xs = sim.get_pvec(0); // this line is optional
		conf.values.ys = plotVecF;
		conf.values.count = sim.get_dim(0);
		
		// run through vector to get min and max in y scale
		float maxValPlot = plotVec[0];
		float minValPlot = plotVec[0];
		for (uint32_t iX = 0; iX < sim.get_dim(0); iX++)
		{
			if (plotVec[iX] > maxValPlot)
				maxValPlot = plotVec[iX];

			if (plotVec[iX] < minValPlot)
				minValPlot = plotVec[iX];
		}

		conf.scale.min = minValPlot;
		conf.scale.max = maxValPlot;
		conf.tooltip.show = true;
		conf.tooltip.format = "x=%.2f, y=%.2f";
		conf.grid_x.show = false;
		conf.grid_y.show = false;
		conf.frame_size = ImVec2(width, 200);
		conf.line_thickness = 2.f;
		ImGui::Plot("plot", conf);

		delete[] plotVecF;

	}


	static char filePath [64];
	ImGui::InputText("Path", filePath, 64);
	string filePathString = filePath;

	ImGui::Columns(2);
	if (ImGui::Button("Export as vtk"))
	{
		bool success = sim.exportVtk(filePathString);
		if (success)
			ImGui::OpenPopup("VTK exported");
		else
			ImGui::OpenPopup("VTK not exported");
      
	}
	ImGui::NextColumn();
	if (ImGui::Button("Export as h5"))
	{
		bool success = sim.exportH5(filePathString);
		if (success)
			ImGui::OpenPopup("H5 exported");
		else
			ImGui::OpenPopup("H5 not exported");      
	}

	ImGui::Columns(1);
	if (ImGui::Button("Export as default  h5"))
	{
		bool success = sim.exportH5();
		if (success)
			ImGui::OpenPopup("H5 exported");
		else
			ImGui::OpenPopup("H5 not exported");      
	}

	if (ImGui::BeginPopupModal("VTK exported"))
  {
    ImGui::Text("Data export successful");
    if (ImGui::Button("Close"))
      ImGui::CloseCurrentPopup();
    ImGui::EndPopup();
  }

  if (ImGui::BeginPopupModal("VTK not exported"))
  {
    ImGui::Text("Data export absolutely not successful");
    if (ImGui::Button("Close"))
      ImGui::CloseCurrentPopup();
    ImGui::EndPopup();
  }

  if (ImGui::BeginPopupModal("H5 exported"))
  {
    ImGui::Text("Data export successful");
    if (ImGui::Button("Close"))
      ImGui::CloseCurrentPopup();
    ImGui::EndPopup();
  }

	if (ImGui::BeginPopupModal("H5 not exported"))
	{
	  ImGui::Text("Data export absolutely not successful");
	  if (ImGui::Button("Close"))
	    ImGui::CloseCurrentPopup();
	  ImGui::EndPopup();
	}


	ImGui::End();
	return;
}


// main loop which we run through to update ImGui
void interface::MainDisplayCode()
{

	Properties(); // TODO: should be removed at some point
	TissueProperties(); // allows defintion of tissue properties
	FieldGenerator(); // allows generation of different fields
	Illumination(); // allows definition of different illumination types
	if (is_output_defined)
	{
		Result();
	}
	return;
}
