#[[
This module defines FetchContent_Declare_Fallback.

FetchContent_Declare_Fallback is a helper function that replaces a call to FetchContent_Declare.
It allows passing options introduced in later CMake versions even if they are not supported
by the current CMake executable.
Currently, the following options are processed:
* FIND_PACKAGE_ARGS (introduced in v3.24)
* SYSTEM (introduced in v3.25)

Note that FIND_PACKAGE_ARGS (and its arguments) has to be passed as the last argument to
FetchContent_Declare_Fallback. This is a limitation of FetchContent_Declare and not specific
to FetchContent_Declare_Fallback.

For more information, consider the CMake documentation for FetchContent_Declare available here:
https://cmake.org/cmake/help/v3.25/module/FetchContent.html#command:fetchcontent_declare
#]]

include_guard(GLOBAL)

include(FetchContent)

function(FetchContent_Declare_Fallback)
    set(FETCH_CONTENT_ARGS ${ARGN})

    if(CMAKE_VERSION VERSION_LESS "3.24.0")
        list(FIND FETCH_CONTENT_ARGS "FIND_PACKAGE_ARGS" FIND_PACKAGE_ARGS_INDEX)

        if(NOT (FIND_PACKAGE_ARGS_INDEX STREQUAL "-1"))
            set(LAST_ELEMENT "")
            while(NOT (LAST_ELEMENT STREQUAL "FIND_PACKAGE_ARGS"))
                list(POP_BACK FETCH_CONTENT_ARGS LAST_ELEMENT)
            endwhile()
        endif()
    endif()

    if(CMAKE_VERSION VERSION_LESS "3.25.0")
        list(REMOVE_ITEM FETCH_CONTENT_ARGS SYSTEM)
    endif()

    FetchContent_Declare(${FETCH_CONTENT_ARGS})
endfunction()
