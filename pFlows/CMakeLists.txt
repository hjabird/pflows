cmake_minimum_required(VERSION 3.1)

# Gather files that constitute library ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
file (GLOB pFlows_INCLUDE "include/pFlows/*.h")
file (GLOB pFlows_SOURCE  "src/*.cpp")
 
source_group("" FILES ${pFlows})
source_group("include" FILES ${pFlows_INCLUDE})
source_group("source" FILES ${pFlows_SOURCE})
 
include_directories (include/pFlows)
 
add_library(pFlows STATIC ${pFlows_INCLUDE} 
                        ${pFlows_SOURCE})

# External library dependencies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# BOOST
find_package(Boost 1.60 REQUIRED)
target_link_libraries(pFlows Boost)
# EIGEN
find_package(Eigen3 REQUIRED)
target_link_libraries(pFlows Eigen)
# HBTK
find_package(HBTK REQUIRED)
message("HBTK " ${HBTK_FOUND})
message("HBTK dir " ${HBTK_DIR})
target_link_libraries(pFlows HBTK)
# target_include_directories(pFlows PRIVATE ${CMAKE_PREFIX_PATH})
 

# Outputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set_property(TARGET pFlows PROPERTY FOLDER "libraries")
 
# Properties->General->Output Directory
set_target_properties(pFlows PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
 
 
# Adds logic to INSTALL.vcproj to copy math.dll to the destination directory
install (TARGETS pFlows
         ARCHIVE DESTINATION ${PROJECT_SOURCE_DIR}/install
         LIBRARY DESTINATION ${PROJECT_SOURCE_DIR}/install)
         
export(PACKAGE pFlows)
