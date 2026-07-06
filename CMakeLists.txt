cmake_minimum_required(VERSION 3.16)
project(particlizer LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Default to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Default the install location to the active conda environment instead of
# /usr/local (which is normally not writable). An explicit
# -DCMAKE_INSTALL_PREFIX=... on the command line always takes precedence.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT AND DEFINED ENV{CONDA_PREFIX})
    set(CMAKE_INSTALL_PREFIX "$ENV{CONDA_PREFIX}" CACHE PATH "Install prefix" FORCE)
    message(STATUS "CMAKE_INSTALL_PREFIX defaulted to active conda env: ${CMAKE_INSTALL_PREFIX}")
endif()

# Compiler optimizations for performance
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -ffast-math -pipe -Wall -Wextra -pedantic")

# Link-time optimization via CMake's IPO machinery instead of a bare -flto in
# CMAKE_CXX_FLAGS. This makes CMake archive static libraries with the compiler's
# plugin-aware gcc-ar/gcc-ranlib, so LTO symbols are written to the archive index
# and the final link resolves them. A plugin-unaware `ar` (e.g. some conda
# binutils) otherwise produces an empty index -> "undefined reference" at link.
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error)
if(ipo_supported)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
else()
    message(WARNING "IPO/LTO not supported, building without it: ${ipo_error}")
endif()

option(BQSSAMPLER_BUILD_EXAMPLES "Build example programs" ON)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# -----------------------------
# Core static library (everything except the CLI main)
# -----------------------------
set(CORE_SOURCES
    src/particle.cpp
    src/particle_system.cpp
    src/surface.cpp
    src/sampler.cpp
    src/numerical_integrator.cpp
    src/tools.cpp
    src/tables.cpp
    src/settings.cpp
)

add_library(sampler_core STATIC ${CORE_SOURCES})
add_library(BQSSampler::sampler_core ALIAS sampler_core)

target_include_directories(sampler_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/bqssampler>
)

target_link_libraries(sampler_core PUBLIC m)

# -----------------------------
# CLI executable (unchanged behavior)
# -----------------------------
add_executable(particlizer src/main.cpp)
target_link_libraries(particlizer PRIVATE sampler_core)

# -----------------------------
# Example / verification driver
# -----------------------------
if(BQSSAMPLER_BUILD_EXAMPLES)
    add_executable(sample_one_cell examples/sample_one_cell.cpp)
    target_link_libraries(sample_one_cell PRIVATE sampler_core)
endif()

# -----------------------------
# Install rules
# -----------------------------
install(TARGETS sampler_core
    EXPORT BQSSamplerTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(TARGETS particlizer
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Public headers
install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/bqssampler
    FILES_MATCHING PATTERN "*.h"
)

# Data tables (assembled grad_coeffs.dat, not the .part_* splits)
install(DIRECTORY tables/
    DESTINATION ${CMAKE_INSTALL_DATADIR}/bqssampler/tables
    FILES_MATCHING PATTERN "*.dat" PATTERN "*.part_*" EXCLUDE
)

# Export targets so consumers can find_package(BQSSampler)
install(EXPORT BQSSamplerTargets
    FILE BQSSamplerTargets.cmake
    NAMESPACE BQSSampler::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/BQSSampler
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/BQSSamplerConfigVersion.cmake"
    VERSION 1.0.0
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/BQSSamplerConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/BQSSamplerConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/BQSSampler
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/BQSSamplerConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/BQSSamplerConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/BQSSampler
)
