
set(OpenGL_GL_PREFERENCE GLVND)
cmake_policy(SET CMP0079 NEW)

# prepare for cuda compilation
enable_language(CUDA)
find_package(CUDA REQUIRED QUIET)
INCLUDE(FindCUDA)
INCLUDE_DIRECTORIES(/opt/cuda/include include)
set(CUDA_SEPARABLE_COMPILATION ON)
set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)
set(BUILD_SHARED_LIBS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# set(CMAKE_CUDA_FLAGS "${CUDA_CUDA_FLAGS} -rdc=true -ftz=false -prec-div=true -prec-sqrt=true")
set(CMAKE_CUDA_FLAGS "${CUDA_CUDA_FLAGS} -rdc=true --ptxas-options=-v -O2")

add_library(Sphere sphere.cpp)
set_property(TARGET Sphere PROPERTY CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(Sphere PUBLIC
	Shape
	)

add_library(Box box.cpp)
set_property(TARGET Box PROPERTY CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(Box PUBLIC
	Shape
	)

add_library(Shell shell.cpp)
set_property(TARGET Shell PROPERTY CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(Shell PUBLIC
	Shape
	)

add_library(Tube tube.cpp)
set_property(TARGET Tube PROPERTY CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(Tube PUBLIC
	Shape
	)

add_library(Shape shape.cpp)
set_property(TARGET Shape PROPERTY CUDA_SEPARABLE_COMPILATION ON)

