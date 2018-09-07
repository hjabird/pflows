cmake_minimum_required(VERSION 3.1)

# Target
add_executable (VortexParticles_demo VortexParticles_demo/VortexParticles_demo.cpp)

# Library dependencies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
target_include_directories (VortexParticles_demo PRIVATE "${PROJECT_SOURCE_DIR}/pFlows/include/pFlows") 
target_link_libraries (VortexParticles_demo pFlows) 
 
# Visual studio ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# VS folders.
set_property(TARGET VortexParticles_demo PROPERTY FOLDER "executables")

# Destinations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set_target_properties (VortexParticles_demo PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# INSTALL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
install (TARGETS VortexParticles_demo
         RUNTIME DESTINATION bin)
