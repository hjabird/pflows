cmake_minimum_required(VERSION 3.1)

# Target
add_executable (Compare2Dand3Dunsteady Compare2Dand3Dunsteady/Compare2Dand3Dunsteady.cpp)

# Library dependencies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
target_include_directories (Compare2Dand3Dunsteady PRIVATE "${PROJECT_SOURCE_DIR}/pFlows/include/pFlows") 
target_link_libraries (Compare2Dand3Dunsteady pFlows)
# BOOST
find_package(Boost 1.60 REQUIRED)
# target_link_libraries(Compare2Dand3Dunsteady Boost)
# EIGEN
# find_package(Eigen3 REQUIRED)
include_directories(Sclavounos1987_demo PUBLIC ${EIGEN3_INCLUDE_DIR})
# target_link_libraries(Compare2Dand3Dunsteady Eigen3) 
 
# Visual studio ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# VS folders.
set_property(TARGET Compare2Dand3Dunsteady PROPERTY FOLDER "executables")

# Destinations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set_target_properties(Compare2Dand3Dunsteady PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# INSTALL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
install (TARGETS Compare2Dand3Dunsteady
         RUNTIME DESTINATION bin)
