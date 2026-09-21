include(CMakePackageConfigHelpers)

set(geosets_config_install_dir "${CMAKE_INSTALL_LIBDIR}/cmake/geosets")

set(geosets_has_vendored_deps FALSE)
foreach(dep_target IN ITEMS libcdd_vendor qhullcpp qhullstatic_r highs eigen boost_geometry)
    if (TARGET ${dep_target})
        get_target_property(is_imported ${dep_target} IMPORTED)
        if (NOT is_imported)
            set(geosets_has_vendored_deps TRUE)
            break()
        endif ()
    endif ()
endforeach()

set(install_export_args)
if (NOT geosets_has_vendored_deps)
    list(APPEND install_export_args EXPORT geosetsTargets)
else ()
    message(STATUS
        "Skipping geosets CMake package export because vendored fallback dependencies are active."
    )
endif ()

install(TARGETS ${cpp_lib_name}
    ${install_export_args}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(
    DIRECTORY ${PROJECT_SOURCE_DIR}/include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.hpp"
    PATTERN ".idea" EXCLUDE
)

if (NOT geosets_has_vendored_deps)
    configure_package_config_file(
        "${PROJECT_SOURCE_DIR}/cmake/geosetsConfig.cmake.in"
        "${PROJECT_BINARY_DIR}/geosetsConfig.cmake"
        INSTALL_DESTINATION "${geosets_config_install_dir}"
    )

    write_basic_package_version_file(
        "${PROJECT_BINARY_DIR}/geosetsConfigVersion.cmake"
        VERSION "${PROJECT_VERSION}"
        COMPATIBILITY SameMajorVersion
    )

    install(
        EXPORT geosetsTargets
        FILE geosetsTargets.cmake
        NAMESPACE geosets::
        DESTINATION "${geosets_config_install_dir}"
    )

    install(
        FILES
        "${PROJECT_BINARY_DIR}/geosetsConfig.cmake"
        "${PROJECT_BINARY_DIR}/geosetsConfigVersion.cmake"
        DESTINATION "${geosets_config_install_dir}"
    )
endif ()
