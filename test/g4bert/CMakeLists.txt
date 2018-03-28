#------------------------------------------------------------------------------
# - CMakeLists.txt for cmsExp
#------------------------------------------------------------------------------
#
# - Define CMake requirements and override make rules as needed
#
cmake_minimum_required(VERSION 2.8)

#------------------------------------------------------------------------------
# - Project definition
#------------------------------------------------------------------------------
project(g4bertTest)

#------------------------------------------------------------------------------
# - Find the package named 'Geant4',
#------------------------------------------------------------------------------
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/include )
FILE(GLOB src_files ${PROJECT_SOURCE_DIR}/src/*.cc )

add_executable(g4bertTest g4bertTest.cc ${src_files})
target_link_libraries(g4bertTest ${Geant4_LIBRARIES})

#------------------------------------------------------------------------------
# End of CMakeLists.txt
#------------------------------------------------------------------------------
