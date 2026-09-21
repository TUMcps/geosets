# Load externalCatch2.cmake script, which finds Catch2
include(externalCatch2)

include(Catch)

# Add all the test files
add_executable(geosets_tests
    definitions.h
    utils.h
    common_tests.cpp
    common_binary_tests.cpp
    parallel_tests.cpp
    sets/interval_tests.cpp
    sets/vpolytope_tests.cpp
    sets/hpolytope_tests.cpp
    sets/polygon_tests.cpp
    sets/zonotope_tests.cpp
    sets/conversion_tests.cpp
    sampling/gaussian_sampling_tests.cpp
)

# Link them with the library
target_link_libraries(geosets_tests
    PRIVATE ${cpp_lib_name} Catch2::Catch2WithMain
)

target_compile_options(geosets_tests PRIVATE
    -Werror=sign-compare
    -Wall
    -Wextra
    -Werror
)

if (GS_USE_OPENMP)
    target_compile_definitions(geosets_tests PRIVATE GS_USE_OPENMP)
endif()

# catch_discover_tests(geosets_tests)
catch_discover_tests(geosets_tests
    PROPERTIES SKIP_RETURN_CODE 4
)

if(GS_BUILD_CODE_COVERAGE)
    message(STATUS "Enabling code coverage instrumentation")

    # Coverage works only with GNU/Clang
    if(NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        message(FATAL_ERROR "Coverage supported only for GCC/Clang.")
    endif()

    # Disable optimizations and add coverage flags
    target_compile_options(${cpp_lib_name} PRIVATE --coverage -O0 -g)
    target_link_options(${cpp_lib_name} PRIVATE --coverage)
    target_compile_options(geosets_tests PRIVATE --coverage -O0 -g)
    target_link_options(geosets_tests PRIVATE --coverage)
endif()
