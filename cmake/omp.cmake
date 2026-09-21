include(FetchContent)

# Workaround to dynamically link to the bundled omp of another python package (e.g. torch)
# This might be fragile
# One could theoretically get around this, by building torch from source with USE_SYSTEM_OPENMP=1, but that takes forever...

# On macOS with Python bindings, prefer torch's bundled libomp if available.
# This avoids the "libomp already initialized" error when both geosets and
# PyTorch are loaded — they must link against the same physical libomp.
set(_omp_lib_path "")
set(_omp_inc_path "")

if (APPLE AND SKBUILD)
    find_package(Python COMPONENTS Interpreter QUIET)
    if (Python_FOUND)
        execute_process(
            COMMAND "${Python_EXECUTABLE}" -c
                "import torch, pathlib; print(pathlib.Path(torch.__file__).parent)"
            OUTPUT_VARIABLE _torch_dir
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
            RESULT_VARIABLE _torch_result
        )
        if (_torch_result EQUAL 0 AND EXISTS "${_torch_dir}/lib/libomp.dylib")
            message(STATUS "Using torch's libomp from: ${_torch_dir}/lib")
            set(_omp_lib_path "${_torch_dir}/lib")
            if (EXISTS "${_torch_dir}/include/omp.h")
                set(_omp_inc_path "${_torch_dir}/include")
            endif ()
        endif ()
    endif ()
endif ()

# Necessary for the AppleClang compiler to find system omp (if installed).
if (APPLE)
    message(STATUS "Looking for OpenMP with apple clang")

    # Fall back to homebrew if torch libomp not found
    if (NOT _omp_lib_path)
        if (${CMAKE_SYSTEM_PROCESSOR} STREQUAL "arm64")
            set(_omp_lib_path /opt/homebrew/opt/libomp/lib)
            set(_omp_inc_path /opt/homebrew/opt/libomp/include)
        else ()
            set(_omp_lib_path /usr/local/opt/libomp/lib)
            set(_omp_inc_path /usr/local/opt/libomp/include)
        endif ()
    elseif (NOT _omp_inc_path)
        # torch libomp found but no headers — use homebrew headers
        if (${CMAKE_SYSTEM_PROCESSOR} STREQUAL "arm64")
            set(_omp_inc_path /opt/homebrew/opt/libomp/include)
        else ()
            set(_omp_inc_path /usr/local/opt/libomp/include)
        endif ()
    endif ()

    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        foreach(_lang IN ITEMS C CXX)
            set(OpenMP_${_lang}_LIB_NAMES "omp")
            set(OpenMP_${_lang}_FLAGS "-Xclang -fopenmp")
            set(OpenMP_${_lang}_INCLUDE_DIR ${_omp_inc_path})
        endforeach()
    endif()

    find_library(OpenMP_omp_LIBRARY
        NAMES omp
        PATHS ${_omp_lib_path}
        NO_DEFAULT_PATH
    )

    set(CMAKE_DISABLE_PRECOMPILE_HEADERS ON)

    # Export the lib path so the python binding can add it as an rpath
    set(GS_OMP_LIB_DIR "${_omp_lib_path}" CACHE INTERNAL "")
endif ()

# Try the system OpenMP first
find_package(OpenMP QUIET)

if (OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found (${OpenMP_CXX_VERSION}): ${OpenMP_CXX_INCLUDE_DIRS}")
else()
    message(STATUS "OpenMP not found. Fetching LLVM OpenMP runtime as fallback...")
    # TODO: it needs to be tested how well this works!

    # Fetch OpenMP runtime source if missing
    FetchContent_Declare(
        llvm_openmp
        GIT_REPOSITORY https://github.com/llvm/llvm-project.git
        GIT_TAG        main
        SOURCE_SUBDIR  openmp
    )

    FetchContent_MakeAvailable(llvm_openmp)

    # Create a pseudo OpenMP::OpenMP_CXX target to maintain API compatibility
    add_library(OpenMP::OpenMP_CXX INTERFACE IMPORTED)
    target_include_directories(OpenMP::OpenMP_CXX INTERFACE
        ${llvm_openmp_SOURCE_DIR}/runtime/src
    )

    # You could add linking rules or compilation flags here
    target_compile_options(OpenMP::OpenMP_CXX INTERFACE -fopenmp)
endif()
