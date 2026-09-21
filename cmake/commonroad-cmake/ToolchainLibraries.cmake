include_guard(GLOBAL)

# Add a proper target for dl, gcov, libatomic, backtrace

# We use LINK_LIBRARIES_ONLY_TARGETS in order to only allow targets to be specified in target_link_libraries.
# However, there is no proper target for libdl and gcov as they is provided by the toolchain (GCC/glibc) itself.
# We create toolchain libraries in order to still have a proper targets for gcov and dl.
# See https://cmake.org/cmake/help/v3.25/prop_tgt/LINK_LIBRARIES_ONLY_TARGETS.html

# dl
add_library(toolchain::dl INTERFACE IMPORTED)

# Note: CMAKE_DL_LIBS may be empty, for example on macOS
if(CMAKE_DL_LIBS)
    set_property(TARGET toolchain::dl PROPERTY IMPORTED_LIBNAME ${CMAKE_DL_LIBS})

    if(NOT TARGET ${CMAKE_DL_LIBS})
        add_library(${CMAKE_DL_LIBS} ALIAS toolchain::dl)
    endif()
endif()

# gcov
add_library(toolchain::gcov INTERFACE IMPORTED)
set_property(TARGET toolchain::gcov PROPERTY IMPORTED_LIBNAME gcov)

if(NOT TARGET gcov)
    add_library(gcov ALIAS toolchain::gcov)
endif()

# libatomic
add_library(toolchain::atomic INTERFACE IMPORTED)

find_library(libatomic atomic)
if(NOT (libatomic MATCHES "libatomic-NOTFOUND"))
    set_property(TARGET toolchain::atomic PROPERTY IMPORTED_LOCATION ${libatomic})
else()
    set_property(TARGET toolchain::atomic PROPERTY IMPORTED_LIBNAME atomic)
endif()

if(NOT TARGET atomic)
    add_library(atomic ALIAS toolchain::atomic)
endif()

# backtrace
find_library(_backtrace_lib backtrace)
add_library(toolchain::backtrace INTERFACE IMPORTED)
if(NOT _backtrace_lib STREQUAL "_backtrace_lib-NOTFOUND")
    message(VERBOSE "found backtrace: ${_backtrace_lib}")
    set_property(TARGET toolchain::backtrace PROPERTY IMPORTED_LOCATION ${_backtrace_lib})
else()
    set_property(TARGET toolchain::backtrace PROPERTY IMPORTED_LIBNAME backtrace)
endif()

if(NOT TARGET backtrace)
    add_library(backtrace ALIAS toolchain::backtrace)
endif()

# rt
add_library(toolchain::rt INTERFACE IMPORTED)
set_property(TARGET toolchain::rt PROPERTY IMPORTED_LIBNAME rt)

if(NOT TARGET rt)
    add_library(rt ALIAS toolchain::rt)
endif()
