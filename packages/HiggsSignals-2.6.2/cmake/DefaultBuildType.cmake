# --------------- set a default build type if none was specified --------------
# defaults to Release
set(default_build_type "Release")
# unless we are in a git checkout of something other than the master
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
  find_package(Git QUIET)
  if(Git_FOUND)
    execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
                    OUTPUT_VARIABLE GIT_HEAD)
    if(NOT ${GIT_HEAD} MATCHES "master")
      set(default_build_type "Debug")
    endif(NOT ${GIT_HEAD} MATCHES "master")
  endif(Git_FOUND)
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(
    STATUS
      "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}"
      CACHE STRING "Choose the type of build."
      FORCE)
  # Set the possible values of build type
  set_property(CACHE CMAKE_BUILD_TYPE
               PROPERTY STRINGS
                        "Debug"
                        "Release"
                        "MinSizeRel"
                        "RelWithDebInfo")
endif()
