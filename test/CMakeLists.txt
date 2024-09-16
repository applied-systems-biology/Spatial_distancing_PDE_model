add_executable(test_configurations  src/testConfigurations.cpp)
target_include_directories(test_configurations PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(test_configurations PRIVATE
        project_options
        abm::simulation
        abm::analyser
        abm::basic
        abm::io
        abm::utils
        Boost::filesystem)

add_executable(test_units  src/testUnits.cpp)
target_include_directories(test_units PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(test_units PRIVATE
        project_options
        abm::analyser
        abm::simulation
        abm::basic
        abm::io
        abm::utils
        Boost::filesystem)


add_test(NAME configurations_functions_tests COMMAND test_configurations)
add_test(NAME analyser_functions_tests COMMAND test_units)