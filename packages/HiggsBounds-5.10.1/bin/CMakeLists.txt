add_executable(HiggsBounds HiggsBounds.F90)
target_link_libraries(HiggsBounds PRIVATE HB)
set_target_properties(HiggsBounds PROPERTIES RUNTIME_OUTPUT_DIRECTORY ..)
add_custom_command(TARGET HiggsBounds POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E create_symlink
                           ${CMAKE_CURRENT_SOURCE_DIR}/../example_data
                           ${CMAKE_CURRENT_BINARY_DIR}/../example_data)

add_executable(AllAnalyses AllAnalyses.F90)
target_link_libraries(AllAnalyses PRIVATE HB)
set_target_properties(AllAnalyses PROPERTIES RUNTIME_OUTPUT_DIRECTORY ..)
