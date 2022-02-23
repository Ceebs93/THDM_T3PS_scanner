find_program(UNAME_PROGRAM NAMES uname)
mark_as_advanced(UNAME_PROGRAM)

if(UNAME_PROGRAM)
  execute_process(COMMAND ${UNAME_PROGRAM} "-m"
                  OUTPUT_VARIABLE MACHINE_NAME
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

if(POLICY CMP0074)
  find_library(FeynHiggs_LIBRARIES FH
               PATH_SUFFIXES "${MACHINE_NAME}-${CMAKE_SYSTEM_NAME}/lib"
                             "${MACHINE_NAME}-${CMAKE_SYSTEM_NAME}/lib64"
                             "build")
  find_path(
    FeynHiggs_INCLUDE_DIRS FHCouplings.h
    PATH_SUFFIXES "${MACHINE_NAME}-${CMAKE_SYSTEM_NAME}/include" "build")
else(POLICY CMP0074)
  find_library(FeynHiggs_LIBRARIES FH
               PATH_SUFFIXES "${MACHINE_NAME}-${CMAKE_SYSTEM_NAME}/lib"
                             "${MACHINE_NAME}-${CMAKE_SYSTEM_NAME}/lib64"
                             "build"
               PATHS ${FeynHiggs_ROOT} ENV FeynHiggs_ROOT)
  find_path(FeynHiggs_INCLUDE_DIRS FHCouplings.h
            PATH_SUFFIXES "${MACHINE_NAME}-${CMAKE_SYSTEM_NAME}/include" "build"
            PATHS ${FeynHiggs_ROOT} ENV FeynHiggs_ROOT)
endif(POLICY CMP0074)

if(FeynHiggs_LIBRARIES AND FeynHiggs_INCLUDE_DIRS)
  if(NOT TARGET FeynHiggs::FH)
    add_library(FeynHiggs::FH UNKNOWN IMPORTED)
    set_target_properties(FeynHiggs::FH
                          PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                     "${FeynHiggs_INCLUDE_DIRS}"
                                     IMPORTED_LOCATION
                                     "${FeynHiggs_LIBRARIES}")
  endif(NOT TARGET FeynHiggs::FH)
endif(FeynHiggs_LIBRARIES AND FeynHiggs_INCLUDE_DIRS)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FeynHiggs
                                  REQUIRED_VARS
                                  FeynHiggs_LIBRARIES
                                  FeynHiggs_INCLUDE_DIRS)
mark_as_advanced(FeynHiggs_LIBRARIES FeynHiggs_INCLUDE_DIRS)
