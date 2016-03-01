cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE  "Release")
INCLUDE_DIRECTORIES(  "$ENV{HOME}/include")
LINK_DIRECTORIES(     "$ENV{HOME}/lib")
SET(CMAKE_CXX_FLAGS   "-Wunused-result -Wno-dev -g -O3")
set(CMAKE_MACOSX_RPATH 1)

#NavierStokes library
project(NavierStokes)
file(                  GLOB_RECURSE source_files sources/*)
file(                  GLOB_RECURSE header_files sources/*.h)
ADD_LIBRARY(           NavierStokes SHARED ${source_files})
find_library(          Commons NAMES Commons HINTS "$ENV{HOME}/lib")
target_link_libraries( NavierStokes ${Commons})
INSTALL(               TARGETS NavierStokes LIBRARY  DESTINATION "$ENV{HOME}/lib")
INSTALL(               FILES ${header_files} DESTINATION "$ENV{HOME}/include")

#Executable
project(main)
add_executable(        nstokes sources/nstokes.c)
target_link_libraries( nstokes Commons NavierStokes)

INSTALL(               TARGETS nstokes RUNTIME DESTINATION "$ENV{HOME}/bin")
set_property(TARGET    nstokes PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)

 



