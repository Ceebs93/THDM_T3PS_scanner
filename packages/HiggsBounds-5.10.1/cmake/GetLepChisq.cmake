find_file(
  LEP_CHISQ_FOUND
  "csboutput_outfile_fullclsb_ee_h2h1_bbbb_090307_1_prodah_dechbb.binary"
  PATHS ${HB_DATA_DIR}/lep-chisq-master/csboutput_trans_binary
  NO_DEFAULT_PATH)

if(NOT LEP_CHISQ_FOUND)
  message(STATUS "LEP chisq tables not found, downloading...")
  file(
    DOWNLOAD
      https://gitlab.com/higgsbounds/lep-chisq/-/archive/master/lep-chisq-master.tar.gz
      ${CMAKE_CURRENT_BINARY_DIR}/lep-chisq.tar.gz
    EXPECTED_HASH MD5=cd6829a80244ae71e9fa4d47ef2c12b3)
  execute_process(COMMAND ${CMAKE_COMMAND} -E tar xfz
                          ${CMAKE_CURRENT_BINARY_DIR}/lep-chisq.tar.gz
                  WORKING_DIRECTORY ${HB_DATA_DIR})
  message(
    STATUS "LEP chisq tables downloaded to ${HB_DATA_DIR}/lep-chisq-master")
endif(NOT LEP_CHISQ_FOUND)

unset(LEP_CHISQ_FOUND CACHE)
