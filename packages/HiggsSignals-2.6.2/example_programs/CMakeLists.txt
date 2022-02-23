file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results)

add_executable(HBandHSwithSLHA HBandHSwithSLHA.f90)
target_link_libraries(HBandHSwithSLHA PRIVATE HS)

add_executable(HS_2Higgses HS_2Higgses.f90)
target_link_libraries(HS_2Higgses PRIVATE HS)

add_executable(HS_efficiencies HS_efficiencies.f90)
target_link_libraries(HS_efficiencies PRIVATE HS)

add_executable(HS_mass HS_mass.f90)
target_link_libraries(HS_mass PRIVATE HS)

add_executable(HS_SM_LHCRun1 HS_SM_LHCRun1.f90)
target_link_libraries(HS_SM_LHCRun1 PRIVATE HS)

add_executable(HSeffC HSeffC.f90)
target_link_libraries(HSeffC PRIVATE HS)

if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU"
   AND "${CMAKE_Fortran_COMPILER_VERSION}" VERSION_GREATER 7.0.0)
  add_executable(HSLimitEffC HSLimitEffC.f90)
  target_link_libraries(HSLimitEffC PRIVATE HS)
endif()

add_executable(HSeffCwithSTXS HSeffCwithSTXS.f90)
target_link_libraries(HSeffCwithSTXS PRIVATE HS)

add_executable(HSgga HSgga.f90)
target_link_libraries(HSgga PRIVATE HS)

add_executable(HShadr HShadr.f90)
target_link_libraries(HShadr PRIVATE HS)

add_executable(HShadrwithSTXS HShadrwithSTXS.f90)
target_link_libraries(HShadrwithSTXS PRIVATE HS)

add_executable(HSwithSLHA HSwithSLHA.f90)
target_link_libraries(HSwithSLHA PRIVATE HS)

add_executable(HSwithSTXS HSwithSTXS.f90)
target_link_libraries(HSwithSTXS PRIVATE HS)

add_executable(HSwithToys HSwithToys.f90)
target_link_libraries(HSwithToys PRIVATE HS)

if(TARGET FeynHiggs::FH)
  add_executable(HBandHSwithFH HBandHSwithFH.F)
  target_link_libraries(HBandHSwithFH PRIVATE HS FeynHiggs::FH)
  target_compile_options(
    HBandHSwithFH
    PRIVATE
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","GNU">:-ffixed-line-length-none>
      $<$<STREQUAL:"${CMAKE_Fortran_COMPILER_ID}","Intel">:-extend_source>)
endif(TARGET FeynHiggs::FH)
