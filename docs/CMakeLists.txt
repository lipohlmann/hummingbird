# Adapted from MFEM, which licensed under BSD-3 https://github.com/mfem/mfem/blob/master/doc/CMakeLists.txt

find_package(Doxygen)

if (DOXYGEN_FOUND)
  set(HUMMINGBIRD_SOURCE_DIR ${PROJECT_SOURCE_DIR})
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

  add_custom_target(doc
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    COMMAND echo "<meta http-equiv=\"REFRESH\" content=\"0;URL=build/html/index.html\">" > ${CMAKE_CURRENT_BINARY_DIR}/documentation.html
    BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/build/html/index.html
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen to documentation.html"
    VERBATIM)

  add_custom_target(clean-doc
    COMMAND ${CMAKE_COMMAND} -E remove -f ${CMAKE_CURRENT_BINARY_DIR}/documentation.html
    COMMAND ${CMAKE_COMMAND} -E remove -f ${CMAKE_CURRENT_BINARY_DIR}/warnings.log
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/build
    COMMENT "Removing API documentation"
    VERBATIM)

endif (DOXYGEN_FOUND)
