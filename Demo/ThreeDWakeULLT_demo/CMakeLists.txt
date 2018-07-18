cmake_minimum_required(VERSION 3.1)

# Target
add_executable (ThreeDWakeULLT_demo ThreeDWakeULLT_demo/ThreeDWakeULLT_demo.cpp)

# Library dependencies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
target_include_directories (ThreeDWakeULLT_demo PRIVATE "${PROJECT_SOURCE_DIR}/pFlows/include/pFlows") 
target_link_libraries (ThreeDWakeULLT_demo pFlows)
# BOOST
find_package(Boost 1.50 REQUIRED)
# EIGEN
find_package(Eigen3 REQUIRED) 
 
# Visual studio ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# VS folders.
set_property(TARGET ThreeDWakeULLT_demo PROPERTY FOLDER "executables")

# Destinations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set_target_properties (ThreeDWakeULLT_demo PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# INSTALL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
install (TARGETS ThreeDWakeULLT_demo
         RUNTIME DESTINATION bin)
