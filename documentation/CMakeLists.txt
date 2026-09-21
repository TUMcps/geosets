find_package(Doxygen REQUIRED)
find_package(Sphinx REQUIRED)

# set input and output files
set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/doxygen/Doxyfile.in)
set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/doxygen/Doxyfile)
configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

set(DOXYGEN_GENERATE_HTML YES)
set(DOXYGEN_GENERATE_XML NO)

# Make sure the python library directory is not included in the C++ documentation
set(DOXYGEN_EXCLUDE_PATTERNS "*/template_py/*")

# For some reason, this does not work in the CI environment
#doxygen_add_docs(doxygen_docs
#    COMMENT "Generate doxygen documentation"
#    CONFIG_FILE ${DOXYGEN_OUT}
#)

# Add Doxygen docs target with config file
add_custom_target(
    doxygen_docs
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generate doxygen documentation"
    DEPENDS "${DOXYGEN_CPPREFERENCE_TAG_FILE}" "${DOXYGEN_OUT}"
    VERBATIM
)

# Add sphinx target with FindSphinx.cmake
sphinx_add_docs(
    sphinx_docs
    BUILDER html
    SOURCE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/sphinx
    OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/sphinx
)
