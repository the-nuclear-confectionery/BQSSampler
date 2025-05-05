cmake_minimum_required(VERSION 3.10)
project(particlizer)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Force Release mode (No Debug mode)
set(CMAKE_BUILD_TYPE Release)

# Compiler optimizations for performance
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -ffast-math -flto -pipe -Wall -Wextra -pedantic")

# Include directories (headers)
include_directories(include)

# Source files
set(SOURCES
    src/main.cpp
    src/particle.cpp
    src/particle_system.cpp
    src/surface.cpp
    src/sampler.cpp
    src/numerical_integrator.cpp
    src/tools.cpp
    src/tables.cpp
    src/settings.cpp
)

# Create the executable
add_executable(particlizer ${SOURCES})

# Link necessary libraries (math library)
target_link_libraries(particlizer m)
