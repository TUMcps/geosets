# Name of the python lib in the .so file
set(python_lib_name "_gs_py_core")

# Find python
set(Python_FIND_VIRTUALENV FIRST)
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
message(STATUS "Python_EXECUTABLE: ${Python_EXECUTABLE}")

# Check if nanobind was already included by a parent project (including it twice leads to a bug)
if (NOT NB_DIR)
    # To circumvent the command being called twice
    execute_process(
        COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
        OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE NB_DIR)
endif ()
list(APPEND CMAKE_PREFIX_PATH "${NB_DIR}")

find_package(nanobind CONFIG REQUIRED)

# Create the nanobind module
nanobind_add_module(${python_lib_name} MODULE
    nanobind_module.cpp
    parallelized_bindings.h
    sampling_bindings.h
)

target_compile_features(${python_lib_name} PRIVATE cxx_std_17)
target_link_libraries(${python_lib_name} PRIVATE ${cpp_lib_name})

if(GS_USE_OPENMP)
    target_compile_definitions(${python_lib_name} PRIVATE GS_USE_OPENMP)

    # Workaround to dynamically link to the bundled omp of another python package (e.g. torch)
    # This might be fragile
    # On macOS, the linked libomp may have an arbitrary install name (e.g.
    # torch ships one with /opt/llvm-openmp/lib/libomp.dylib). Rewrite it
    # to @rpath/libomp.dylib and add the actual directory as an rpath, so
    # the .so finds whichever libomp is available at runtime.
    if(APPLE AND GS_OMP_LIB_DIR)
        set_target_properties(${python_lib_name} PROPERTIES
            BUILD_RPATH "${GS_OMP_LIB_DIR}"
            INSTALL_RPATH "${GS_OMP_LIB_DIR}"
        )
        # Get the install name of the libomp we linked against
        execute_process(
            COMMAND otool -D "${GS_OMP_LIB_DIR}/libomp.dylib"
            OUTPUT_VARIABLE _omp_install_name
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        # otool -D output: first line is the file path, second is the install name
        string(REGEX MATCH "\n(.*)" _match "${_omp_install_name}")
        set(_omp_id "${CMAKE_MATCH_1}")
        if(_omp_id AND NOT _omp_id STREQUAL "@rpath/libomp.dylib")
            add_custom_command(TARGET ${python_lib_name} POST_BUILD
                COMMAND install_name_tool -change "${_omp_id}" "@rpath/libomp.dylib"
                        "$<TARGET_FILE:${python_lib_name}>"
                COMMENT "Rewriting libomp install name for macOS compatibility"
            )
        endif()
    endif()
endif()

nanobind_add_stub(
    geosets_stub
    MODULE ${python_lib_name}
    OUTPUT ${python_lib_name}.pyi
    MARKER_FILE py.typed
    PYTHON_PATH $<TARGET_FILE_DIR:${python_lib_name}>
    DEPENDS ${python_lib_name}
    RECURSIVE true
    INCLUDE_PRIVATE true  # Necessary to include functions with ending _
)

# Install the library into the python environment folder
if (SKBUILD)
    message(STATUS "The python bindings are built using scikit-build-core")
    install(TARGETS ${python_lib_name}
        LIBRARY DESTINATION ${SKBUILD_PROJECT_NAME}
    )

    # New layout py.typed + _gs_py_core/__init__.pyi + _gs_py_core/sampling.pyi
    # Install stubs
    install(
        FILES ${CMAKE_CURRENT_BINARY_DIR}/py.typed
        DESTINATION ${SKBUILD_PROJECT_NAME}
    )

    # Install the stub package directory
    install(
        DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${python_lib_name}
        DESTINATION ${SKBUILD_PROJECT_NAME}
    )

    # Install into src directory as well (for better development)
    install(
        FILES ${CMAKE_CURRENT_BINARY_DIR}/py.typed
        DESTINATION ${PROJECT_SOURCE_DIR}/src/${SKBUILD_PROJECT_NAME}
    )
    install(
        DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${python_lib_name}
        DESTINATION ${PROJECT_SOURCE_DIR}/src/${SKBUILD_PROJECT_NAME}
    )
endif ()
