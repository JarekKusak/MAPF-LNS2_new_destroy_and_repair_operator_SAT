cmake_minimum_required(VERSION 3.16)
project(lns CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# Optional: suppress warnings (pozdeji vymazat)
set(CMAKE_CXX_FLAGS "-w")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Zdroje a hlavičky
file(GLOB_RECURSE SRC
        src/*.cpp
        src/CBS/*.cpp
        src/PIBT/*.cpp)

add_executable(lns ${SRC})

target_include_directories(lns PUBLIC
        inc inc/CBS inc/PIBT
        include
        include/libs/monosat
        include/libs/monosat/utils
        include/libs/monosat/graph)

#  Third‑party balíčky
find_package(Boost 1.70 REQUIRED COMPONENTS program_options system filesystem)
target_link_libraries(lns PRIVATE
        Boost::program_options Boost::system Boost::filesystem)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(lns PRIVATE Eigen3::Eigen)

if(APPLE)
    # Homebrew fallback, když by find_package selhal
    find_library(GMP_LIB   NAMES gmp   PATHS /opt/homebrew/Cellar/gmp/*/lib)
    find_library(GMPXX_LIB NAMES gmpxx PATHS /opt/homebrew/Cellar/gmp/*/lib)
    find_library(ZLIB_LIB  NAMES z     PATHS /opt/homebrew/Cellar/zlib/*/lib)
    target_link_libraries(lns PRIVATE ${GMP_LIB} ${GMPXX_LIB} ${ZLIB_LIB})
elseif(UNIX)
    message(STATUS "Detected Linux – using system libraries for GMP and Zlib")

    # ====== najdeme hlavičky a statické/ sdílené knihovny přímo ======
    find_path(GMP_INCLUDE_DIR  gmp.h)
    find_library(GMP_LIB       NAMES gmp)
    find_library(GMPXX_LIB     NAMES gmpxx)
    find_library(ZLIB_LIB      NAMES z)

    if(NOT GMP_INCLUDE_DIR OR NOT GMP_LIB OR NOT GMPXX_LIB OR NOT ZLIB_LIB)
        message(FATAL_ERROR "GMP / Zlib development files nejsou nainstalovány")
    endif()

    include_directories(${GMP_INCLUDE_DIR})
    set(GMP_LIBRARIES  ${GMP_LIB}  ${GMPXX_LIB})
    set(ZLIB_LIBRARIES ${ZLIB_LIB})
endif()

# Statické knihovny z MAPF‑encodings
foreach(lib IN ITEMS mapf pb kissat monosat)
    add_library(${lib}_lib STATIC IMPORTED)
    set_target_properties(${lib}_lib
            PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/libs/lib${lib}.a)
endforeach()

# Link‑order: nejdřív knihovna, která symboly _potřebuje_ (mapf),
# až pak její závislosti.
target_link_libraries(lns PRIVATE mapf_lib pb_lib kissat_lib monosat_lib)