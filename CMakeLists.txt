cmake_minimum_required(VERSION 3.10)
project(particlizer)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Enable debug symbols for debugging (optional)
set(CMAKE_BUILD_TYPE Debug)

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
)

# Add executable
add_executable(particlizer ${SOURCES})

# Link necessary libraries (if applicable)
# Example: link math library
target_link_libraries(particlizer m)
