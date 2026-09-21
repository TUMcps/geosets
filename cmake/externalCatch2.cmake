include(FetchContent)

# Trying to find Catch2 on the system
find_package(Catch2 QUIET)
if(Catch2_FOUND)
    message(STATUS "Catch2 found (${Catch2_VERSION}): ${Catch2_INCLUDE_DIRS}")
else ()
    message(STATUS "Catch2 not found, fetching fallback (${Catch2_VERSION}) ...")
    FetchContent_Declare(Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.4.0 # or a later release
    )
    FetchContent_MakeAvailable(Catch2)
endif()
