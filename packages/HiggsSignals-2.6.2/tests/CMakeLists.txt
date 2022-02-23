list(
  APPEND
  MHMODPLUS_FILES
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
    COMMAND
      ${CMAKE_COMMAND} -E create_symlink
      "${CMAKE_CURRENT_SOURCE_DIR}/../example_data/mhmodplus/${file}"
      "${CMAKE_CURRENT_BINARY_DIR}/mhmodplus_effC/${file}")
  execute_process(
    COMMAND
      ${CMAKE_COMMAND} -E create_symlink
      "${CMAKE_CURRENT_SOURCE_DIR}/../example_data/mhmodplus/${file}"
      "${CMAKE_CURRENT_BINARY_DIR}/mhmodplus_hadr/${file}")
endforeach(file IN LISTS MHMODPLUS_FILES)

add_test(
  NAME "HiggsSignals_effC"
  COMMAND
    ${PYTHON_EXECUTABLE} -m pytest
    ${CMAKE_CURRENT_SOURCE_DIR}/test_HiggsSignals.py::test_HiggsSignals_effC
  WORKING_DIRECTORY $<TARGET_FILE_DIR:HiggsSignals>)

add_test(
  NAME "HiggsSignals_hadr"
  COMMAND
    ${PYTHON_EXECUTABLE} -m pytest
    ${CMAKE_CURRENT_SOURCE_DIR}/test_HiggsSignals.py::test_HiggsSignals_hadr
  WORKING_DIRECTORY $<TARGET_FILE_DIR:HiggsSignals>)

add_test(
  NAME "HS_2Higgses"
  COMMAND ${PYTHON_EXECUTABLE} -m pytest
          ${CMAKE_CURRENT_SOURCE_DIR}/test_HiggsSignals.py::test_HS_2Higgses
  WORKING_DIRECTORY $<TARGET_FILE_DIR:HS_2Higgses>)

add_test(
  NAME "HSeffC"
  COMMAND ${PYTHON_EXECUTABLE} -m pytest
          ${CMAKE_CURRENT_SOURCE_DIR}/test_HiggsSignals.py::test_HSeffC
  WORKING_DIRECTORY $<TARGET_FILE_DIR:HSeffC>)

add_test(
  NAME "HShadr"
  COMMAND ${PYTHON_EXECUTABLE} -m pytest
          ${CMAKE_CURRENT_SOURCE_DIR}/test_HiggsSignals.py::test_HShadr
  WORKING_DIRECTORY $<TARGET_FILE_DIR:HShadr>)

file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/../example_data/SLHA/SLHA_FHexample.fh.1"
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(
  NAME "HSwithSLHA"
  COMMAND ${PYTHON_EXECUTABLE} -m pytest
          ${CMAKE_CURRENT_SOURCE_DIR}/test_HiggsSignals.py::test_HSwithSLHA
  WORKING_DIRECTORY $<TARGET_FILE_DIR:HSwithSLHA>)

add_test(
  NAME "HS_mass"
  COMMAND ${PYTHON_EXECUTABLE} -m pytest
          ${CMAKE_CURRENT_SOURCE_DIR}/test_HiggsSignals.py::test_HS_mass
  WORKING_DIRECTORY $<TARGET_FILE_DIR:HS_mass>)

add_test(
  NAME "HS_SM_LHCRun1"
  COMMAND ${PYTHON_EXECUTABLE} -m pytest
          ${CMAKE_CURRENT_SOURCE_DIR}/test_HiggsSignals.py::test_HS_SM_LHCRun1
  WORKING_DIRECTORY $<TARGET_FILE_DIR:HS_SM_LHCRun1>)

if(TARGET HSLimitEffC)
  add_test(
    NAME "HSLimitEffC"
    COMMAND ${PYTHON_EXECUTABLE} -m pytest
            ${CMAKE_CURRENT_SOURCE_DIR}/test_HiggsSignals.py::test_HSLimitEffC
    WORKING_DIRECTORY $<TARGET_FILE_DIR:HSLimitEffC>)
endif()
