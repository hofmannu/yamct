# yamct
Here we go, yet another toolbox for light simulations in highly scattering media like biological tissues. Completely unnecessary one might think. There are so many out there someone else might mention. Nevertheless, when I started to look for such a toolbox I was highly disappointed by the non-userfriendly versions which people published in the web over the last few years. Either using them is impossible due to a bad user interface or they perform incredibly slow because they are still running on a CPU.

This is my attempt to build a Monte Carlo simulation toolbox featuring CUDA acceleration on NVIDIA cards, a friendly user interface based on the amazing `ImGui` project, handling multilayered tissues with elegance, geometry creation on the fly and directly from the GUI, and export to datatypes readable from different software like Python, MATLAB, etc. through the established `h5` standard. The simulations are performed voxel by voxel relying on Code extensively tested on the CPU.

![Main interface of YAMCT](https://hofmannu.org/wp-content/uploads/2020/12/yamct-1024x795.jpg "Main interface")

# Installation

The code was written fro and tested with Linux. Please send me a request if you do not know how to install it on your system and I will provide you with instructions.

## Installation of required libraries on Linux

This part is still fully missing since I did not have the time yet to check which libraries I installed over the last half a year. Things we definetly need are
*  `hdf5` general purpose library and file format for storing scientific data
*  `cuda` and `nvidia` driver (eventually you want to use `nvidia-lts` if you use `linux-lts`)
*   stuff required for imgui
*  `cmake` and `make`

Archlinux installation command:
`pacman -S cuda nvidia hdf5 cmake make`

## Building and running the program
Please do not continue here if you are not using Linux. This software is built for Linux machines and requires furthermore CUDA. Switch on the terminal, `cd` to your favorite installation directory and run the following command cascade:

```
git clone git@github.com:hofmannu/yamct.git
cd yamct
git submodule init
git submodule update
mkdir Debug
cd Debug
cmake .. && make all && ./main_exp
```

# Support

I am actively working on this project. If you want any feature implemented (for example different geometrical shapes, export types, or illumination types) feel free to open an issue and I will get back to you as soon as possible.

# Further reading

# Similar / alternative projects
