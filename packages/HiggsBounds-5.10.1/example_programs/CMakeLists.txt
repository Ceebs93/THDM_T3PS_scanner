add_executable(HBeffC HBeffC.F90)
target_link_libraries(HBeffC PRIVATE HB)

add_executable(HBwithchannelrates HBwithchannelrates.F90)
target_link_libraries(HBwithchannelrates PRIVATE HB)

if(TARGET FeynHiggs::FH)
  add_executable(HBwithFH HBwithFH.F)
  target_link_libraries(HBwithFH PRIVATE HB FeynHiggs::FH)

  add_executable(HBSLHAinputblocksfromFH HBSLHAinputblocksfromFH.F90
                 HBSLHAinputblocksfromFH_extras.F)
  target_link_libraries(HBSLHAinputblocksfromFH PRIVATE HB FeynHiggs::FH)

  foreach(TARGET_WITH_FH IN ITEMS HBwithFH HBSLHAinputblocksfromFH)
    target_compile_options(
      ${TARGET_WITH_FH}
      PRIVATE
        $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-ffixed-line-length-none>
        $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:-extend_source>)
  endforeach()
endif(TARGET FeynHiggs::FH)

if(LEP_CHISQ)
  add_executable(HBwithLEPlikelihood HBwithLEPlikelihood.F90)
  target_link_libraries(HBwithLEPlikelihood PRIVATE HB)
endif(LEP_CHISQ)

add_executable(HBwithLHClikelihood HBwithLHClikelihood.F90)
target_link_libraries(HBwithLHClikelihood PRIVATE HB)

add_executable(HBwithSLHA HBwithSLHA.F90)
target_link_libraries(HBwithSLHA PRIVATE HB)
