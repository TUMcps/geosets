#### Project setup ####
cmake_minimum_required(VERSION 3.16)
set(CMAKE_CXX_STANDARD 17)

# Setup CMake metadata
set(project_name "geosets")
set(project_version "0.1.0")

# Retrieve metadata from pyproject.toml, if the project is built with scikit-build
if (SKBUILD)
    set(project_name ${SKBUILD_PROJECT_NAME})
    set(project_version ${SKBUILD_PROJECT_VERSION})
endif ()

project(
    ${project_name}
    VERSION ${project_version}
    DESCRIPTION "C++ geometric sets with Python bindings."
    LANGUAGES C CXX
)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
include(GNUInstallDirs)

# Add cmake folder to the module path
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# Define the name of the cpp library
set(cpp_lib_name "geosets_lib")

#### Options ####
option(GS_BUILD_CPP_EXAMPLE "Build the cpp example" OFF)
option(GS_BUILD_PYTHON_BINDINGS "Build the python bindings" OFF)
option(GS_BUILD_TESTS "Build the cpp tests" OFF)
option(GS_BUILD_DOCUMENTATION "Build the cpp documentation with doxygen" OFF)
option(GS_USE_OPENMP "Use OpenMP for parallelized computations" OFF)

include(CMakeDependentOption) # For dependent options
cmake_dependent_option(
    GS_BUILD_CODE_COVERAGE "Build code coverage" OFF
    "GS_BUILD_TESTS;NOT SKBUILD" OFF
)

message(STATUS "GS_BUILD_CPP_EXAMPLE: ${GS_BUILD_CPP_EXAMPLE}")
message(STATUS "GS_BUILD_PYTHON_BINDINGS: ${GS_BUILD_PYTHON_BINDINGS}")
message(STATUS "GS_BUILD_TESTS: ${GS_BUILD_TESTS}")
message(STATUS "GS_BUILD_DOCUMENTATION: ${GS_BUILD_DOCUMENTATION}")
message(STATUS "GS_BUILD_CODE_COVERAGE: ${GS_BUILD_CODE_COVERAGE}")
message(STATUS "GS_USE_OPENMP: ${GS_USE_OPENMP}")

add_subdirectory(src)

if (GS_BUILD_CPP_EXAMPLE)
    add_subdirectory(examples/cpp)
endif ()

if (GS_BUILD_PYTHON_BINDINGS)
    add_subdirectory(bindings/python)
endif ()

if (GS_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests/cpp)
endif ()

if (GS_BUILD_DOCUMENTATION)
    add_subdirectory(documentation)
endif ()

include(geosetsPackage)
