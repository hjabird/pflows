cmake_minimum_required(VERSION 3.1)

# Target
add_executable (Ramesh2014_demo Ramesh2014_demo/Ramesh2014_demo.cpp)

# Library dependencies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
target_include_directories (Ramesh2014_demo PRIVATE "${PROJECT_SOURCE_DIR}/pFlows/include/pFlows") 
target_link_libraries (Ramesh2014_demo pFlows)
# BOOST
find_package(Boost 1.60 REQUIRED)
# target_link_libraries(Sclavounos1987_demo Boost)
# EIGEN
find_package(Eigen3 REQUIRED)
# target_link_libraries(Sclavounos1987_demo Eigen3) 
 
# Visual studio ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# VS folders.
set_property(TARGET Ramesh2014_demo PROPERTY FOLDER "executables")

# Destinations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set_target_properties(Ramesh2014_demo PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# INSTALL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
install (TARGETS Ramesh2014_demo
         RUNTIME DESTINATION bin)
