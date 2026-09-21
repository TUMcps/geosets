# Library files

set(header_files
    ${PROJECT_SOURCE_DIR}/include/geosets/core.h
    ${PROJECT_SOURCE_DIR}/include/geosets/structs.h
    ${PROJECT_SOURCE_DIR}/include/geosets/linalg_utils.h
    ${PROJECT_SOURCE_DIR}/include/geosets/exceptions.h
    ${PROJECT_SOURCE_DIR}/include/geosets/parallel.h
    ${PROJECT_SOURCE_DIR}/include/geosets/qhull_wrapper.h
    ${PROJECT_SOURCE_DIR}/include/geosets/highs_wrapper.h
    ${PROJECT_SOURCE_DIR}/include/geosets/infinite.h
    ${PROJECT_SOURCE_DIR}/include/geosets/sampling/common.h
    ${PROJECT_SOURCE_DIR}/include/geosets/sampling/gaussian_combined.h
    ${PROJECT_SOURCE_DIR}/include/geosets/sampling/gaussian_rdhr.h
    ${PROJECT_SOURCE_DIR}/include/geosets/sampling/gaussian_rejection.h

    ${PROJECT_SOURCE_DIR}/include/geosets/conversion/to_hpolytope.h
    ${PROJECT_SOURCE_DIR}/include/geosets/conversion/to_polygon.h
    ${PROJECT_SOURCE_DIR}/include/geosets/conversion/to_vpolytope.h

    ${PROJECT_SOURCE_DIR}/include/geosets/sets/interface.h
    ${PROJECT_SOURCE_DIR}/include/geosets/sets/hpolytope.h
    ${PROJECT_SOURCE_DIR}/include/geosets/sets/vpolytope.h
    ${PROJECT_SOURCE_DIR}/include/geosets/sets/polygon.h
    ${PROJECT_SOURCE_DIR}/include/geosets/sets/interval.h
    ${PROJECT_SOURCE_DIR}/include/geosets/sets/zonotope.h
)

set(source_files
    linalg_utils.cpp
    highs_wrapper.cpp
    sampling/gaussian_common.cpp
    sampling/gaussian_combined.cpp
    sampling/gaussian_rdhr.cpp
    sampling/gaussian_rejection.cpp

    sets/hpolytope.cpp
    sets/vpolytope.cpp
    sets/polygon.cpp
    sets/interval.cpp
    sets/zonotope.cpp

    conversion/to_hpolytope.cpp
    conversion/to_polygon.cpp
    conversion/to_vpolytope.cpp
)

# External packages as defined from commonroad cmake (removed, since the repo is private)
# include(commonroadCMake)

# External package include is defined in cmake/eigen3.cmake
include(eigen3)
include(boost)
include(highs)
include(qhull)
include(libcdd)

# Create the library
add_library(${cpp_lib_name}
    ${source_files} ${header_files}
)
set_target_properties(${cpp_lib_name} PROPERTIES EXPORT_NAME geosets)
add_library(geosets::geosets ALIAS ${cpp_lib_name})


target_compile_features(${cpp_lib_name} PRIVATE cxx_std_17)

target_link_libraries(${cpp_lib_name}
    PUBLIC
    Eigen3::Eigen
    Boost::geometry
    PRIVATE
    Qhull::qhullcpp
    highs::highs
    libcdd::cdd
)

# More strict error handling
target_compile_options(${cpp_lib_name} PRIVATE
    -Wall
    -Wextra
    -Werror
)

# Optional include and linking of openmp
if (GS_USE_OPENMP)
    include(omp)
    target_link_libraries(${cpp_lib_name} PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(${cpp_lib_name} PRIVATE GS_USE_OPENMP)
endif()

target_include_directories(${cpp_lib_name}
    PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
