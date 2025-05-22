cmake_minimum_required(VERSION 3.16)
project(lns CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Výchozí build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

# Volitelně zapnout AddressSanitizer
option(ENABLE_ASAN "Enable AddressSanitizer" OFF)
if (ENABLE_ASAN)
    message(STATUS "AddressSanitizer enabled")
    add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
    add_link_options(-fsanitize=address -fno-omit-frame-pointer)
endif()

# Soubory projektu
file(GLOB_RECURSE SRC
    src/*.cpp
    src/CBS/*.cpp
    src/PIBT/*.cpp)

add_executable(lns ${SRC})

# Hlavičky
target_include_directories(lns PUBLIC
    inc inc/CBS inc/PIBT
    include
    include/libs/monosat
    include/libs/monosat/utils
    include/libs/monosat/graph)

# Boost
find_package(Boost 1.70 REQUIRED COMPONENTS program_options system filesystem)
target_link_libraries(lns PRIVATE Boost::program_options Boost::system Boost::filesystem)

# Eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(lns PRIVATE Eigen3::Eigen)

# GMP a Zlib – Linux vs macOS
if(APPLE)
    find_library(GMP_LIB     NAMES gmp   PATHS /opt/homebrew/Cellar/gmp/*/lib REQUIRED)
    find_library(GMPXX_LIB   NAMES gmpxx PATHS /opt/homebrew/Cellar/gmp/*/lib REQUIRED)
    find_library(ZLIB_LIB    NAMES z     PATHS /opt/homebrew/Cellar/zlib/*/lib REQUIRED)
    find_path  (GMP_INCLUDE  NAMES gmp.h PATHS /opt/homebrew/Cellar/gmp/*/include REQUIRED)
else()
    find_library(GMP_LIB     NAMES gmp   REQUIRED)
    find_library(GMPXX_LIB   NAMES gmpxx REQUIRED)
    find_path  (GMP_INCLUDE  NAMES gmp.h REQUIRED)
    find_package(ZLIB REQUIRED)
endif()

target_include_directories(lns PRIVATE ${GMP_INCLUDE})

# Statické knihovny z MAPF-encodings
foreach(lib IN ITEMS mapf pb kissat monosat)
    if(NOT TARGET ${lib}_lib)
        add_library(${lib}_lib STATIC IMPORTED)
        set_target_properties(${lib}_lib PROPERTIES IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/libs/lib${lib}.a)
        target_link_libraries(lns PRIVATE ${lib}_lib)
    endif()
endforeach()

# Linknutí systémových knihoven až PO statických
target_link_libraries(lns PRIVATE ${GMP_LIB} ${GMPXX_LIB} ZLIB::ZLIB)
