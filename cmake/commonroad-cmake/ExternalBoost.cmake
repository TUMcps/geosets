include(FetchContent)

# CMP0167 - FindBoost removed
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()

FetchContent_Declare(
  Boost
  URL https://github.com/boostorg/boost/releases/download/boost-1.87.0/boost-1.87.0-cmake.tar.xz
  URL_HASH SHA256=7da75f171837577a52bbf217e17f8ea576c7c246e4594d617bfde7fafd408be5
  DOWNLOAD_EXTRACT_TIMESTAMP ON
  EXCLUDE_FROM_ALL
  SYSTEM

  FIND_PACKAGE_ARGS 1.87 COMPONENTS ${BOOST_INCLUDE_LIBRARIES} CONFIG

  # OPTIONAL_COMPONENTS program_options filesystem CONFIG
)

# Specify the libraries we need in order to reduce compile times
# NOTE: The Boost docs mention that you should determine transitive dependencies
# by using `boostdep --brief`, however this does not appear to be necesary
#
# This should be defined by the consumers of this module, however you need to ensure that
# the set of libraries you include here covers all requirements of the project itself and its dependencies
# (as Boost can only be imported once)
#
# set(BOOST_INCLUDE_LIBRARIES geometry program_options filesystem stacktrace)

# For dl, backtrace, libatomic
include(commonroad-cmake/ToolchainLibraries)

FetchContent_MakeAvailable(Boost)

# Add more header-only targets here as required
foreach(target geometry graph align polygon)
    if(NOT TARGET Boost::${target})
        add_library(Boost::${target} ALIAS Boost::headers)
    endif()
endforeach()

# Fix Boost::stacktrace CMake configuration
# (they export an include directory that references the build directory)
foreach(target boost_stacktrace_backtrace boost_stacktrace_addr2line)
    if(TARGET ${target})
        # get_target_property(includes ${target} INTERFACE_INCLUDE_DIRECTORIES)
        # message(STATUS "Boost includes: ${includes}")
        set_property(TARGET ${target} PROPERTY INTERFACE_INCLUDE_DIRECTORIES
          $<BUILD_INTERFACE:${includes}>
          $<INSTALL_INTERFACE:include>
        )
    endif()
endforeach()

if(TARGET boost_stacktrace_backtrace)
    set_property(TARGET boost_stacktrace_backtrace PROPERTY LINK_LIBRARIES_ONLY_TARGETS OFF)
endif()
