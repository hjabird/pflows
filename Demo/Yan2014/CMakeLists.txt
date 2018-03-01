cmake_minimum_required(VERSION 3.1)

# Target
add_executable (Yan2014 Yan2014/Yan2014.cpp)

# Library dependencies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
target_include_directories (Yan2014 PRIVATE "${PROJECT_SOURCE_DIR}/pFlows/include/pFlows") 
target_link_libraries (Yan2014 pFlows)
# BOOST
find_package(Boost 1.60 REQUIRED)
# target_link_libraries(Yan2014 Boost)
# EIGEN
find_package(Eigen3 REQUIRED)
# target_link_libraries(Yan2014 Eigen3) 
 
# Visual studio ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# VS folders.
set_property(TARGET Yan2014 PROPERTY FOLDER "executables")

# Destinations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set_target_properties(Yan2014 PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# INSTALL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
install (TARGETS Yan2014
         RUNTIME DESTINATION bin)
