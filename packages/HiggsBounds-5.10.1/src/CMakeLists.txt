add_library(
  HB
  ${CMAKE_CURRENT_BINARY_DIR}/install_data.f90
  usefulbits.f90
  S95tables_type1.F90
  S95tables_type2.F90
  S95tables_type3.F90
  interpolate.f90
  theory_colliderSfunctions.f90
  theory_XS_SM_functions.f90
  theory_BRfunctions.F90
  likelihoods.F90
  S95tables.f90
  PDGnumbering.f90
  string_manip.f90
  SLHA_manip.f90
  extra_bits_for_SLHA.f90
  input.F90
  channels.F90
  output.F90
  theo_manip.f90
  HiggsBounds_subroutines.F90
  access_SM.f90
  access_effC.f90
  CInterface/access_effC_C.f90
  CInterface/access_SM_C.f90
  CInterface/HiggsBounds_subroutines_C.F90)
set_target_properties(HB PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ../lib
                                    LIBRARY_OUTPUT_DIRECTORY ../lib)
target_compile_options(
  HB PRIVATE $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-Wall -Wextra -Wno-surprising>
             $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:-warn>)
target_include_directories(HB PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
                                     ${PROJECT_SOURCE_DIR}/include)
if(LEP_CHISQ)
  target_compile_definitions(HB PRIVATE enableCHISQ)
  target_sources(HB PRIVATE extra_bits_for_chisquared.F90)
endif(LEP_CHISQ)
add_custom_command(
  TARGET HB POST_BUILD COMMAND ${CMAKE_COMMAND} -DHB_DATA_DIR=${HB_DATA_DIR} -P
                               ${PROJECT_SOURCE_DIR}/cmake/RefreshBindata.cmake)
