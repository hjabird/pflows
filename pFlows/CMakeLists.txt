cmake_minimum_required(VERSION 3.1)

# Collect sources into the variable SIMPLE_FUNCTION_SOURCES
file (GLOB pFlows_INCLUDE "include/pFlows/*.h")
file (GLOB pFlows_SOURCE  "src/*.cpp")
 
# Create named folders for the sources within the .vcproj
# Empty name lists them directly under the .vcproj
source_group("" FILES ${pFlows})
source_group("include" FILES ${pFlows_INCLUDE})
source_group("source" FILES ${pFlows_SOURCE})
 
# Properties->C/C++->General->Additional Include Directories
include_directories (include/pFlows)
 
# Set Properties->General->Configuration Type to Dynamic Library(.dll)
# Creates math.dll with the listed sources collected in the variables
# Also adds sources to the Solution Explorer
add_library(pFlows STATIC ${pFlows_INCLUDE} 
                        ${pFlows_SOURCE})
find_package(HBTK)
target_link_libraries(pFlows HBTK::HBTK)
 
# Creates folder "libraries" and adds target project (math.vcproj)
set_property(TARGET pFlows PROPERTY FOLDER "libraries")
 
# Properties->General->Output Directory
set_target_properties(pFlows PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
 
 
# Adds logic to INSTALL.vcproj to copy math.dll to the destination directory
install (TARGETS pFlows
         ARCHIVE DESTINATION ${PROJECT_SOURCE_DIR}/install
         LIBRARY DESTINATION ${PROJECT_SOURCE_DIR}/install)
         
export(PACKAGE pFlows)
