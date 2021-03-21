/*  
	interfacing class to use our volumetric monte carlo simulation

	Author: Urs Hofmann
	Mail: hofmannu@ethz.ch
	Date: 04.05.2020

	ToDo
		- allow multiple illumination types at the same time
		- implement automatic loading for tissue types from files
		- allow automatic loading of illumination types from files

	Changelog

*/

#ifndef INTERFACE_H
#define INTERFACE_H

#include <SDL2/SDL.h>
#include <GL/glew.h>    // Initialize with gl3wInit()

#include <thread> 

#include "../lib/imgui/imgui.h"
#include "imgui_impl_sdl.h"
#include "imgui_impl_opengl3.h"
#include "ImGuiFileDialog.h"
#include "imgui_plot.h"
#include "mc.cuh"
#include "fiberProperties.h"
#include "mcFieldProperties.h"
#include "optProperties.h"
#include "simProperties.h"
#include <cstdio>
#include "color_mapper.cuh"

// different geometrical shapes
#include "optVolume.h"
#include "sphere.h"
#include "box.h"	
#include "tube.h"

using namespace std;

class interface
{
public:
	void InitWindow(int *argcp, char**argv);
	interface();
	~interface();

private:
	void ScanGpus(); // function updates internal GPU list
	void MainDisplayCode();
	void Properties();
	void Result(); // display result to user
	void TissueProperties(); // define tissue properties
	void FieldGenerator(); // define volume of absorbers and scatterers
	void Illumination();
	void ImImagesc(
	const float* data, const uint64_t sizex, const uint64_t sizey, 
	GLuint* out_texture, const color_mapper myCMap);

	void ImImagesc(
	const double* data, const uint64_t sizex, const uint64_t sizey, 
	GLuint* out_texture, const color_mapper myCMap);

	const char* windowTitle = "YAMCT";
	// ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 0.10f);
	ImVec4 clear_color = ImVec4(0.60f, 0.55f, 0.45f, 0.10f);


	bool show_properties_window = 1;
	bool show_tissue_properties = 1;
	bool show_results = 1;
	bool show_field_generator = 1;
	bool show_illumination = 1;

	bool is_output_defined = 0;
	bool is_volume_generated = 0;
	bool is_simulated = 0;

	// slice position in x.y.z for output volume
	float xPos = 0;
	float yPos = 0;
	float zPos = 0;

	bool flagLog = 1;
	bool flagFluence = 1;

	int deviceCount = 0;
	vector<cudaDeviceProp> deviceNames;

	mc sim;
	vector<fiberProperties>* fibers;
	simProperties* simprop;

	vector<optProperties>* tissueTypes; // vector containing different tissue types
	optVolume* volume; // pointer will be stolen from sim

	// elements required for plotting of resulting fluence
	color_mapper fluence_mapper;
	GLuint fluence_texture_x;
	GLuint fluence_texture_y;

	color_mapper material_mapper;
	GLuint material_texture;

	bool flagLogPlot = 1;

};

#endif