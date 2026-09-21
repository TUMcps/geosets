include(FetchContent)

set(_LIBCDD_READY FALSE)

# 1) Config-package installs
find_package(cddlib CONFIG QUIET)
if (TARGET cddlib::cdd)
    add_library(libcdd::cdd ALIAS cddlib::cdd)
    message(STATUS "libcdd found: cddlib::cdd")
    set(_LIBCDD_READY TRUE)
elseif (TARGET cdd)
    add_library(libcdd::cdd ALIAS cdd)
    message(STATUS "libcdd found: cdd")
    set(_LIBCDD_READY TRUE)
endif ()

# 2) Traditional find module installs
if (NOT _LIBCDD_READY)
    find_package(CDD QUIET)
    if (TARGET CDD::CDD)
        add_library(libcdd::cdd ALIAS CDD::CDD)
        message(STATUS "libcdd found: CDD::CDD")
        set(_LIBCDD_READY TRUE)
    elseif (CDD_FOUND)
        add_library(libcdd::cdd INTERFACE IMPORTED)
        if (DEFINED CDD_INCLUDE_DIRS)
            set_target_properties(libcdd::cdd PROPERTIES
                INTERFACE_INCLUDE_DIRECTORIES "${CDD_INCLUDE_DIRS}"
            )
        endif ()
        if (DEFINED CDD_LIBRARIES)
            set_target_properties(libcdd::cdd PROPERTIES
                INTERFACE_LINK_LIBRARIES "${CDD_LIBRARIES}"
            )
        elseif (DEFINED CDD_LIBRARY)
            set_target_properties(libcdd::cdd PROPERTIES
                INTERFACE_LINK_LIBRARIES "${CDD_LIBRARY}"
            )
        endif ()
        message(STATUS "libcdd found via CDD variables")
        set(_LIBCDD_READY TRUE)
    endif ()
endif ()

# 3) pkg-config installs
if (NOT _LIBCDD_READY)
    find_package(PkgConfig QUIET)
    if (PkgConfig_FOUND)
        pkg_check_modules(CDDLIB QUIET cddlib)
        if (CDDLIB_FOUND)
            add_library(libcdd::cdd INTERFACE IMPORTED)
            set_target_properties(libcdd::cdd PROPERTIES
                INTERFACE_INCLUDE_DIRECTORIES "${CDDLIB_INCLUDE_DIRS}"
                INTERFACE_LINK_LIBRARIES "${CDDLIB_LINK_LIBRARIES}"
            )
            message(STATUS "libcdd found via pkg-config")
            set(_LIBCDD_READY TRUE)
        endif ()
    endif ()
endif ()

# 4) find_library/find_path fallback (e.g. Homebrew without pkg-config)
if (NOT _LIBCDD_READY)
    find_library(CDD_LIBRARY_PATH NAMES cdd
        HINTS /opt/homebrew/lib /usr/local/lib)
    find_path(CDD_INCLUDE_PATH NAMES cddlib/cdd.h
        HINTS /opt/homebrew/include /usr/local/include)
    if (CDD_LIBRARY_PATH AND CDD_INCLUDE_PATH)
        add_library(libcdd::cdd INTERFACE IMPORTED)
        set_target_properties(libcdd::cdd PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${CDD_INCLUDE_PATH}"
            INTERFACE_LINK_LIBRARIES "${CDD_LIBRARY_PATH}"
        )
        message(STATUS "libcdd found via find_library: ${CDD_LIBRARY_PATH}")
        set(_LIBCDD_READY TRUE)
    endif ()
endif ()

# 5) Fetch source and build minimal target in-tree
if (NOT _LIBCDD_READY)
    message(STATUS "libcdd not found, fetching fallback...")
    FetchContent_Declare(
        libcdd
        GIT_REPOSITORY https://github.com/cddlib/cddlib.git
        GIT_TAG master
    )
    FetchContent_MakeAvailable(libcdd)

    add_library(libcdd_vendor STATIC
        ${libcdd_SOURCE_DIR}/lib-src/cddcore.c
        ${libcdd_SOURCE_DIR}/lib-src/cddio.c
        ${libcdd_SOURCE_DIR}/lib-src/cddlib.c
        ${libcdd_SOURCE_DIR}/lib-src/cddlp.c
        ${libcdd_SOURCE_DIR}/lib-src/cddmp.c
        ${libcdd_SOURCE_DIR}/lib-src/cddproj.c
        ${libcdd_SOURCE_DIR}/lib-src/setoper.c
    )
    # Copy headers into a cddlib/ subdir so #include <cddlib/cdd.h> works
    set(_LIBCDD_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_libcdd_include")
    file(COPY "${libcdd_SOURCE_DIR}/lib-src/"
        DESTINATION "${_LIBCDD_INCLUDE_DIR}/cddlib"
        FILES_MATCHING PATTERN "*.h"
    )
    target_include_directories(libcdd_vendor PUBLIC "${_LIBCDD_INCLUDE_DIR}")
    # Suppress warnings from third-party libcdd source code
    target_compile_options(libcdd_vendor PRIVATE -w)
    add_library(libcdd::cdd ALIAS libcdd_vendor)
    message(STATUS "libcdd fallback target created from fetched source")
endif ()
