cmake_minimum_required(VERSION 3.16)
project(lns CXX)

# Tool‑chain / language settings
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-w") # optional: silence warnings

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# LIBS directory depends on OS:  libs-macos‑macos/  vs  libs-macos‑linux/
if(APPLE)
    set(LNS_LIBDIR "${CMAKE_SOURCE_DIR}/libs-macos")
else()
    set(LNS_LIBDIR "${CMAKE_SOURCE_DIR}/libs-linux")
endif()

#  Project sources / headers
file(GLOB_RECURSE SRC
     src/*.cpp
     src/CBS/*.cpp
     src/PIBT/*.cpp)

add_executable(lns ${SRC})

target_include_directories(lns PUBLIC
    inc          inc/CBS     inc/PIBT
    include
    include/libs/monosat
    include/libs/monosat/utils
    include/libs/monosat/graph)

# Third‑party (Boost, Eigen, GMP, Zlib)
find_package(Boost 1.70 REQUIRED COMPONENTS program_options system filesystem)
target_link_libraries(lns PRIVATE Boost::program_options Boost::system Boost::filesystem)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(lns PRIVATE Eigen3::Eigen)

if(APPLE)
    # Homebrew fallback (when find_package(GMP) fails on macOS)
    find_library(GMP_LIB     NAMES gmp   PATHS /opt/homebrew/Cellar/gmp/*/lib     REQUIRED)
    find_library(GMPXX_LIB   NAMES gmpxx PATHS /opt/homebrew/Cellar/gmp/*/lib     REQUIRED)
    find_library(ZLIB_LIB    NAMES z     PATHS /opt/homebrew/Cellar/zlib/*/lib    REQUIRED)
    find_path   (GMP_INCLUDE NAMES gmp.h PATHS /opt/homebrew/Cellar/gmp/*/include REQUIRED)
else()
    find_library(GMP_LIB     NAMES gmp     REQUIRED)
    find_library(GMPXX_LIB   NAMES gmpxx   REQUIRED)
    find_path   (GMP_INCLUDE NAMES gmp.h   REQUIRED)
    find_package(ZLIB REQUIRED) # provides ZLIB::ZLIB
endif()

target_include_directories(lns PRIVATE ${GMP_INCLUDE})

# Static libraries produced by MAPF‑encodings
# (libmapf.a, libpb.a, libkissat.a, libmonosat.a)
foreach(lib IN ITEMS mapf pb kissat monosat)
    if(NOT TARGET ${lib}_lib)
        add_library(${lib}_lib STATIC IMPORTED)
        set_target_properties(${lib}_lib PROPERTIES
            IMPORTED_LOCATION "${LNS_LIBDIR}/lib${lib}.a")
        # Static libs-macos first → they *use* the symbols from the system libs-macos
        target_link_libraries(lns PRIVATE ${lib}_lib)
    endif()
endforeach()

# Link system libs-macos **after** the static ones
if(APPLE)
    target_link_libraries(lns PRIVATE ${GMP_LIB} ${GMPXX_LIB} ${ZLIB_LIB})
else()
    target_link_libraries(lns PRIVATE ${GMP_LIB} ${GMPXX_LIB} ZLIB::ZLIB)
endif()