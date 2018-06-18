cmake_minimum_required(VERSION 3.1)

# Target
add_executable (PrandtlLiftingLine_demo PrandtlLiftingLine_demo/PrandtlLiftingLine_demo.cpp)

# Library dependencies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
target_include_directories (PrandtlLiftingLine_demo PRIVATE "${PROJECT_SOURCE_DIR}/pFlows/include/pFlows") 
target_link_libraries (PrandtlLiftingLine_demo pFlows)
# BOOST
find_package(Boost 1.50 REQUIRED)
# EIGEN
find_package(Eigen3 REQUIRED) 
 
# Visual studio ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# VS folders.
set_property(TARGET PrandtlLiftingLine_demo PROPERTY FOLDER "executables")

# Destinations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set_target_properties (PrandtlLiftingLine_demo PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# INSTALL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
install (TARGETS PrandtlLiftingLine_demo
         RUNTIME DESTINATION bin)
