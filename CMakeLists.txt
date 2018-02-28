cmake_minimum_required(VERSION 3.1)
project(pFlows)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
 
# Turn on the ability to create folders to organize projects (.vcproj)
# It creates "CMakePredefinedTargets" folder by default and adds CMake
# defined projects like INSTALL.vcproj and ZERO_CHECK.vcproj
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
 
# Set compiler flags and options. 
# Here it is setting the Visual Studio warning level to 3
# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3")
# Command to output information to the console
# Useful for displaying errors, warnings, and debugging
message ("cxx Flags:" ${CMAKE_CXX_FLAGS})
 
# Sub-directories where more CMakeLists.txt exist
# add_subdirectory(Demo)
add_subdirectory (pFlows)
 
# Turn on CMake testing capabilities
# enable_testing()
 
# Add test cases
#add_test(AppTest1 ${PROJECT_BINARY_DIR}/bin/app 100)
#add_test(AppTest2 ${PROJECT_BINARY_DIR}/bin/app 200)
#add_test(AppTest3 ${PROJECT_BINARY_DIR}/bin/app 300)