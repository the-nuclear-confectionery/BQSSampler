# deltaf_coeff_generator/CMakeLists.txt
cmake_minimum_required(VERSION 3.16)
project(deltaf_coeff_generator LANGUAGES CXX)

# ---------- C++ standard ----------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---------- Paths ----------
# This CMakeLists lives in deltaf_coeff_generator/
set(PROJECT_ROOT "${CMAKE_CURRENT_LIST_DIR}/..")
set(INC_DIR       "${PROJECT_ROOT}/include")
set(SRC_DIR       "${PROJECT_ROOT}/src")

# ---------- Executable ----------
add_executable(build_df_coeff
  "${CMAKE_CURRENT_LIST_DIR}/main.cpp"

  # Core sources used by your generator main:
  "${SRC_DIR}/numerical_integrator.cpp"
  "${SRC_DIR}/particle_system.cpp"
  "${SRC_DIR}/settings.cpp"
  "${SRC_DIR}/tables.cpp"
  "${SRC_DIR}/tools.cpp"

  # Add more here ONLY if the linker complains you need them:
  # "${SRC_DIR}/surface.cpp"
  # "${SRC_DIR}/sampler.cpp"
)

target_include_directories(build_df_coeff PRIVATE
  "${INC_DIR}"
  "${CMAKE_CURRENT_LIST_DIR}"   # so it can see deltaf_integrals.h etc if local
)

# ---------- Warnings (optional but useful) ----------
target_compile_options(build_df_coeff PRIVATE
  -O3
  -march=native
  -Wall -Wextra -Wpedantic
)

# If you use OpenMP later, uncomment:
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(build_df_coeff PRIVATE OpenMP::OpenMP_CXX)
endif()

# If you ever need libm explicitly (usually not needed with g++)
# target_link_libraries(build_df_coeff PRIVATE m)
