cmake_minimum_required(VERSION 3.3)
project(
  HiggsSignals
  VERSION 2.6.2
  LANGUAGES C Fortran)

# -------------------------------- setup cmake --------------------------------
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW) # use <packagename>_ROOT as search path
endif(POLICY CMP0074)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(DefaultBuildType)

# ----------------------- create install_data_HS module -----------------------
set(HS_DATA_DIR
    "${CMAKE_CURRENT_SOURCE_DIR}/data"
    CACHE PATH "HiggsSignals data directory")
string(LENGTH ${HS_DATA_DIR} path_length)
configure_file(src/install_data_HS.f90.in src/install_data_HS.f90 @ONLY)

# -------------------------------- dependencies -------------------------------
if(NOT TARGET HiggsBounds::HB)
  find_package(HiggsBounds 5.10.0 REQUIRED)
endif()
if(PROJECT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  find_package(FeynHiggs)
endif()

# ----------------------------------- doc -------------------------------------
if(PROJECT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  include(MakeDoxygen)
endif()

# ------------------------------ subdirectories -------------------------------
add_subdirectory(src)
if(PROJECT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  add_subdirectory(example_programs)
endif()

# ----------------------- alias for use as subdirectory -----------------------
add_library(HiggsSignals::HS ALIAS HS)

# -------------------------- export as cmake package --------------------------
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  HiggsSignalsConfigVersion.cmake
  VERSION ${PACKAGE_VERSION}
  COMPATIBILITY SameMajorVersion)

export(
  TARGETS HS
  NAMESPACE HiggsSignals::
  FILE HiggsSignalsConfig.cmake)
if(PROJECT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  export(PACKAGE HiggsBounds)
endif()
# ----------------------------------- ctest -----------------------------------
if(PROJECT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  find_package(PythonInterp 3.5.0)
  if(${PythonInterp_FOUND})
    enable_testing()
    add_subdirectory("tests")
  endif()
endif()
