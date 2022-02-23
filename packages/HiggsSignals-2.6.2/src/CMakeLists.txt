add_library(
  HS
  ${CMAKE_CURRENT_BINARY_DIR}/install_data_HS.f90
  all_chisq.f90
  CInterface/HiggsSignals_subroutines_C.f90
  combinatorics.f90
  datatables.f90
  evaluate.f90
  expt_syst.f90
  HiggsSignals_subroutines.F90
  io.f90
  mc_chisq.f90
  numerics.f90
  pc_chisq.f90
  STXS.f90
  usefulbits_HS.F90)
target_link_libraries(HS PUBLIC HiggsBounds::HB)
target_compile_options(
  HS PRIVATE $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-Wall -Wextra>
             $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:-warn>)
set_target_properties(HS PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ../lib
                                    LIBRARY_OUTPUT_DIRECTORY ../lib)
target_include_directories(HS PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
                                     ${PROJECT_SOURCE_DIR}/include)

add_executable(HiggsSignals HiggsSignals.F90)
target_link_libraries(HiggsSignals PRIVATE HS)
set_target_properties(HiggsSignals PROPERTIES RUNTIME_OUTPUT_DIRECTORY ..)
add_custom_command(
  TARGET HiggsSignals POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${CMAKE_CURRENT_SOURCE_DIR}/../example_data
          ${CMAKE_CURRENT_BINARY_DIR}/../example_data)
