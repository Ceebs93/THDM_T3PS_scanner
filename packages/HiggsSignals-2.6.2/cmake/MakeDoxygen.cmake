find_package(Doxygen)

if(DOXYGEN_FOUND)
  if(DOXYGEN_MATHJAX_RELPATH)
    set(MATHJAX_RELPATH "MATHJAX_RELPATH = ${DOXYGEN_MATHJAX_RELPATH}")
  endif()
  configure_file(${PROJECT_SOURCE_DIR}/doc/Doxyfile.in Doxyfile @ONLY)
  add_custom_target(
    doc
    COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generate Documentation with Doxygen")
endif()
