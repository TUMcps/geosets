include(FetchContent)
find_package(Qhull QUIET)
if (Qhull_FOUND)
    message(STATUS "Qhull found (${Qhull_VERSION}): ${Qhull_INCLUDE_DIRS}")

    # qhullcpp depends on the C reentrant library but upstream's CMake config
    # doesn't declare the transitive dependency when qhullcpp is static.
    get_target_property(_qhullcpp_type Qhull::qhullcpp TYPE)
    if (_qhullcpp_type STREQUAL "STATIC_LIBRARY")
        if (BUILD_SHARED_LIBS)
            set(_qhull_c Qhull::qhull_r)
        else()
            set(_qhull_c Qhull::qhullstatic_r)
        endif()
        set_property(TARGET Qhull::qhullcpp APPEND PROPERTY
            INTERFACE_LINK_LIBRARIES ${_qhull_c})
    endif()
else ()
    message(STATUS "Qhull not found, fetching fallback (${Qhull_VERSION}) ...")
    FetchContent_Declare(
        qhull
        GIT_REPOSITORY https://github.com/qhull/qhull.git
        GIT_TAG        master
    )
    FetchContent_MakeAvailable(qhull)

    if (NOT TARGET Qhull::qhullcpp)
        add_library(Qhull::qhullcpp ALIAS qhullcpp)
    endif()
    target_link_libraries(qhullcpp PUBLIC qhullstatic_r)
endif ()
