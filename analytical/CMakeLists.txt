cmake_minimum_required(VERSION 3.16)
project(analytical_cf LANGUAGES CXX)

# -----------------------------
# C++ standard
# -----------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# -----------------------------
# Default build type
# -----------------------------
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# -----------------------------
# OpenMP (for CPU parallelization)
# -----------------------------
find_package(OpenMP REQUIRED)

# -----------------------------
# Project structure
# -----------------------------
set(PROJ_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/..)
set(INC_DIR    ${PROJ_ROOT}/include)
set(SRC_DIR    ${PROJ_ROOT}/src)

# -----------------------------
# Core sources
# -----------------------------
set(CORE_SOURCES
  ${SRC_DIR}/surface.cpp
  ${SRC_DIR}/particle_system.cpp
  ${SRC_DIR}/particle.cpp
  ${SRC_DIR}/numerical_integrator.cpp
  ${SRC_DIR}/settings.cpp
  ${SRC_DIR}/tools.cpp
  ${SRC_DIR}/tables.cpp
  ${SRC_DIR}/sampler.cpp
)

# -----------------------------
# Core static library
# -----------------------------
add_library(sampler_core STATIC ${CORE_SOURCES})

target_include_directories(sampler_core PUBLIC ${INC_DIR})

target_compile_options(sampler_core PRIVATE
  -O3
  -march=native
  -Wall
  -Wextra
  -Wpedantic
)

# Link OpenMP here so everything using sampler_core inherits it
target_link_libraries(sampler_core PUBLIC OpenMP::OpenMP_CXX)

# -----------------------------
# Analytical executable
# -----------------------------
add_executable(ana main_ana.cpp)

target_include_directories(ana PRIVATE ${INC_DIR})

target_compile_options(ana PRIVATE
  -O3
  -march=native
)

target_link_libraries(ana PRIVATE sampler_core)

# -----------------------------
# Clean binary location
# -----------------------------
set_target_properties(ana PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)
