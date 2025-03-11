project(lns)

cmake_minimum_required(VERSION 3.16)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

# Include directories
include_directories(
        "inc" "inc/CBS" "inc/PIBT" "include"
        "include/libs/monosat" "include/libs/monosat/utils" "include/libs/monosat/graph"
        /opt/homebrew/Cellar/gmp/6.3.0/include
        /opt/homebrew/Cellar/zlib/1.3.1/include
)

# Add source files
file(GLOB SOURCES "src/*.cpp" "src/CBS/*.cpp" "src/PIBT/*.cpp")
add_executable(lns ${SOURCES} src/LNS.cpp)

# Add path to libraries
link_directories(${CMAKE_SOURCE_DIR}/include/libs)

# Link static libraries from MAPF
set(MAPF_LIBRARIES
        ${CMAKE_SOURCE_DIR}/libs/libpb.a
        ${CMAKE_SOURCE_DIR}/libs/libkissat.a
        ${CMAKE_SOURCE_DIR}/libs/libmonosat.a
        ${CMAKE_SOURCE_DIR}/libs/libmapf.a
)

# Find Boost
find_package(Boost REQUIRED COMPONENTS program_options system filesystem)

# Find Eigen3 for PIBT
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Add linker flags for GMP and Zlib
set(GMP_LIBRARIES /opt/homebrew/Cellar/gmp/6.3.0/lib/libgmp.a /opt/homebrew/Cellar/gmp/6.3.0/lib/libgmpxx.a)
set(ZLIB_LIBRARIES /opt/homebrew/Cellar/zlib/1.3.1/lib/libz.a)
set(CMAKE_CXX_FLAGS "-w") # později vymazat

# Link all libraries
target_link_libraries(lns ${Boost_LIBRARIES} Eigen3::Eigen ${GMP_LIBRARIES} ${ZLIB_LIBRARIES} ${MAPF_LIBRARIES})

# Print build settings for verification
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Boost Libraries: ${Boost_LIBRARIES}")
message(STATUS "Eigen3 Include Dir: ${Eigen3_INCLUDE_DIRS}")
message(STATUS "GMP Libraries: ${GMP_LIBRARIES}")
message(STATUS "Zlib Libraries: ${ZLIB_LIBRARIES}")
message(STATUS "MAPF Libraries: ${MAPF_LIBRARIES}")
