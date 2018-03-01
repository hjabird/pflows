cmake_minimum_required(VERSION 3.1)

# Target
add_executable (Sclavounos1987_demo Sclavounos1987/Sclavounos1987.cpp)

# Library dependencies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
target_include_directories (Sclavounos1987_demo PRIVATE "${PROJECT_SOURCE_DIR}/pFlows/include/pFlows") 
target_link_libraries (Sclavounos1987_demo pFlows)
# BOOST
find_package(Boost 1.60 REQUIRED)
# target_link_libraries(Sclavounos1987_demo Boost)
# EIGEN
find_package(Eigen3 REQUIRED)
# target_link_libraries(Sclavounos1987_demo Eigen3) 
 
# Visual studio ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# VS folders.
set_property(TARGET Sclavounos1987_demo PROPERTY FOLDER "executables")

# Destinations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set_target_properties(Sclavounos1987_demo PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# INSTALL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
install (TARGETS Sclavounos1987_demo
         RUNTIME DESTINATION bin)
