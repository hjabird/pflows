cmake_minimum_required(VERSION 3.1)

# Target
add_executable (PlanarWakeULLT_demo PlanarWakeULLT_demo/PlanarWakeULLT_demo.cpp)

# Library dependencies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
target_include_directories (PlanarWakeULLT_demo PRIVATE "${PROJECT_SOURCE_DIR}/pFlows/include/pFlows") 
target_link_libraries (PlanarWakeULLT_demo pFlows)
# BOOST
find_package(Boost 1.60 REQUIRED)
# EIGEN
find_package(Eigen3 REQUIRED) 
 
# Visual studio ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# VS folders.
set_property(TARGET PlanarWakeULLT_demo PROPERTY FOLDER "executables")

# Destinations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set_target_properties (PlanarWakeULLT_demo PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# INSTALL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
install (TARGETS PlanarWakeULLT_demo
         RUNTIME DESTINATION bin)
