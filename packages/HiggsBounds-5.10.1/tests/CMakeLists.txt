add_executable(CInterfaceTest CInterfaceTest.c)
target_link_libraries(CInterfaceTest PRIVATE HB)

list(
  APPEND MHMODPLUS_FILES
         "mhmod+_additional.dat"
         "mhmod+_BR_H_NP.dat"
         "mhmod+_BR_H_OP.dat"
         "mhmod+_BR_Hplus.dat"
         "mhmod+_BR_t.dat"
         "mhmod+_CP_values.dat"
         "mhmod+_effC.dat"
         "mhmod+_LEP_2H_CS_ratios.dat"
         "mhmod+_LEP_H_ff_CS_ratios.dat"
         "mhmod+_LEP_HpHm_CS_ratios.dat"
         "mhmod+_LEP_HZ_CS_ratios.dat"
         "mhmod+_LHC13_1H_hadCS_ratios.dat"
         "mhmod+_LHC13_Hplus_hadCS.dat"
         "mhmod+_LHC7_1H_hadCS_ratios.dat"
         "mhmod+_LHC8_1H_hadCS_ratios.dat"
         "mhmod+_LHC8_Hplus_hadCS.dat"
         "mhmod+_MHall_uncertainties.dat"
         "mhmod+_MH_GammaTot.dat"
         "mhmod+_MHplus_GammaTot.dat"
         "mhmod+_TEV_1H_hadCS_ratios.dat")

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mhmodplus_effC")
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mhmodplus_hadr")

foreach(file IN LISTS MHMODPLUS_FILES)
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E create_symlink
            "${CMAKE_CURRENT_SOURCE_DIR}/../example_data/mhmodplus/${file}"
            "${CMAKE_CURRENT_BINARY_DIR}/mhmodplus_effC/${file}")
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E create_symlink
            "${CMAKE_CURRENT_SOURCE_DIR}/../example_data/mhmodplus/${file}"
            "${CMAKE_CURRENT_BINARY_DIR}/mhmodplus_hadr/${file}")
endforeach(file IN LISTS MHMODPLUS_FILES)

add_test(
  NAME "HiggsBounds_effC"
  COMMAND ${PYTHON_EXECUTABLE} -m pytest
          ${CMAKE_CURRENT_SOURCE_DIR}/test_HiggsBounds.py::test_HiggsBounds_effC
  WORKING_DIRECTORY $<TARGET_FILE_DIR:HiggsBounds>)
add_test(
  NAME "HiggsBounds_hadr"
  COMMAND ${PYTHON_EXECUTABLE} -m pytest
          ${CMAKE_CURRENT_SOURCE_DIR}/test_HiggsBounds.py::test_HiggsBounds_hadr
  WORKING_DIRECTORY $<TARGET_FILE_DIR:HiggsBounds>)
add_test(
  NAME "HBeffC"
  COMMAND ${PYTHON_EXECUTABLE} -m pytest
          ${CMAKE_CURRENT_SOURCE_DIR}/test_HiggsBounds.py::test_HBeffC
  WORKING_DIRECTORY $<TARGET_FILE_DIR:HBeffC>)
add_test(
  NAME "HBwithchannelrates"
  COMMAND
    ${PYTHON_EXECUTABLE} -m pytest
    ${CMAKE_CURRENT_SOURCE_DIR}/test_HiggsBounds.py::test_HBwithchannelrates
  WORKING_DIRECTORY $<TARGET_FILE_DIR:HBwithchannelrates>)
add_test(
  NAME "HBwithLHClikelihood"
  COMMAND
    ${PYTHON_EXECUTABLE} -m pytest
    ${CMAKE_CURRENT_SOURCE_DIR}/test_HiggsBounds.py::test_HBwithLHClikelihood
  WORKING_DIRECTORY $<TARGET_FILE_DIR:HBwithLHClikelihood>)

if(LEP_CHISQ)
  add_test(
    NAME "HBwithLEPlikelihood"
    COMMAND
      ${PYTHON_EXECUTABLE} -m pytest
      ${CMAKE_CURRENT_SOURCE_DIR}/test_HiggsBounds.py::test_HBwithLEPlikelihood
    WORKING_DIRECTORY $<TARGET_FILE_DIR:HBwithLEPlikelihood>)
endif(LEP_CHISQ)

add_executable(InvisibleBR T_InvisibleBR.f90)
target_link_libraries(InvisibleBR PRIVATE HB)
add_test(InvisibleBR ${CMAKE_CURRENT_BINARY_DIR}/InvisibleBR)
