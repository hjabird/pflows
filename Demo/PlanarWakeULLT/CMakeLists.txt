cmake_minimum_required(VERSION 3.1)

# Target
add_executable (PlanarWakeULLT PlanarWakeULLT/PlanarWakeULLT.cpp)

# Library dependencies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
target_include_directories (PlanarWakeULLT PRIVATE "${PROJECT_SOURCE_DIR}/pFlows/include/pFlows") 
target_link_libraries (PlanarWakeULLT pFlows)
# BOOST
find_package(Boost 1.60 REQUIRED)
# target_link_libraries(Sclavounos1987_demo Boost)
# EIGEN
find_package(Eigen3 REQUIRED)
# target_link_libraries(Sclavounos1987_demo Eigen3) 
 
# Visual studio ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# VS folders.
set_property(TARGET PlanarWakeULLT PROPERTY FOLDER "executables")

# Destinations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set_target_properties(PlanarWakeULLT PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# INSTALL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
install (TARGETS PlanarWakeULLT
         RUNTIME DESTINATION bin)
