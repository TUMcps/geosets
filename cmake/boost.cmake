include(FetchContent)

if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()

set(BOOST_INCLUDE_LIBRARIES geometry)
set(boost_ver 1.87.0)

FetchContent_Declare(
    Boost
    URL https://github.com/boostorg/boost/releases/download/boost-${boost_ver}/boost-${boost_ver}-cmake.tar.xz
    URL_HASH SHA256=7da75f171837577a52bbf217e17f8ea576c7c246e4594d617bfde7fafd408be5
    DOWNLOAD_EXTRACT_TIMESTAMP ON
    EXCLUDE_FROM_ALL
    SYSTEM
    FIND_PACKAGE_ARGS 1.81 CONFIG
)
FetchContent_MakeAvailable(Boost)

if(NOT TARGET Boost::geometry)
    add_library(Boost::geometry ALIAS Boost::headers)
endif()

if(Boost_FOUND)
    message(STATUS "Boost ${Boost_VERSION} found: ${Boost_INCLUDE_DIRS}")
else()
    message(STATUS "Boost ${boost_ver} fetched as fallback")
endif()
