find_package(HiGHS QUIET)

if(HiGHS_FOUND)
    # Upstream HiGHS leaks bare "/include" into INTERFACE_INCLUDE_DIRECTORIES
    get_target_property(_HIGHS_INCS highs::highs INTERFACE_INCLUDE_DIRECTORIES)
    list(REMOVE_ITEM _HIGHS_INCS "/include")
    set_target_properties(highs::highs PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${_HIGHS_INCS}"
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_HIGHS_INCS}"
    )
    message(STATUS "HiGHS found (${HiGHS_VERSION}). include dirs: ${_HIGHS_INCS}")
else ()
    message(STATUS "HiGHS not found, fetching fallback (${HiGHS_VERSION}) ...")
    FetchContent_Declare(
        highs
        GIT_REPOSITORY https://github.com/ERGO-Code/HiGHS.git
        GIT_TAG master
    )

    # Build options for faster & lighter integration
    set(HIGHS_BUILD_SHARED OFF CACHE BOOL "" FORCE)
    set(HIGHS_BUILD_CLI OFF CACHE BOOL "" FORCE)
    set(HIGHS_BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(HIGHS_DEV OFF CACHE BOOL "" FORCE)

    # Necessary to override the shared lib option when building python bindings
    if (SKBUILD)
        set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
    endif ()

    FetchContent_MakeAvailable(highs)
    set(HiGHS_FOUND TRUE)

    set_target_properties(highs PROPERTIES
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
        "$<TARGET_PROPERTY:highs,INTERFACE_INCLUDE_DIRECTORIES>"
    )
endif()
