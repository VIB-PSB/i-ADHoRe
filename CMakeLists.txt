# CMake compatibility issues: don't modify this, please!
cmake_minimum_required(VERSION 2.6)

# user section
project(i-ADHoRe)

# project version
set(${PROJECT_NAME}_MAJOR_VERSION 2)
set(${PROJECT_NAME}_MINOR_VERSION 4)
set(${PROJECT_NAME}_PATCH_LEVEL 1)

# set the module path
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# build static libraries
option(BUILD_SHARED_LIBS OFF)

# set default cmake build type to Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

# make sure we have threading support
find_package(Threads REQUIRED)

# configure MPI
set(MPI_LIBRARY "")
set(MPI_LIBRARIES "")
find_package(MPI)
if (MPI_FOUND)
    add_definitions(-DHAVE_MPI)
    include_directories(${MPI_INCLUDE_PATH})
    set(CMAKE_SHARED_LINKER_FLAGS ${MPI_LINK_FLAGS})
endif(MPI_FOUND)

# configure PNG
find_package(PNG)
if (PNG_FOUND)
    add_definitions(-DHAVE_PNG)
    include_directories(${PNG_INCLUDE_DIR})
    add_definitions(${PNG_DEFINITIONS})
endif(PNG_FOUND)

# set the target-specific flags
set(CMAKE_CXX_FLAGS "-Wno-deprecated ${MPI_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g3")

# set(CMAKE_VERBOSE_MAKEFILE ON)

# configure install
install(DIRECTORY API DESTINATION .)
install(DIRECTORY post_processing DESTINATION .)
install(DIRECTORY testset DESTINATION .)

add_subdirectory(src)
add_subdirectory(unittest)
