add_subdirectory(analyser)
add_subdirectory(basic)
add_subdirectory(io)
add_subdirectory(simulation)
add_subdirectory(utils)

target_link_libraries(analyser PRIVATE project_options project_warnings OpenMP::OpenMP_CXX)
target_link_libraries(basic PRIVATE project_options project_warnings OpenMP::OpenMP_CXX)
target_link_libraries(io PRIVATE project_options project_warnings OpenMP::OpenMP_CXX)
target_link_libraries(simulation PRIVATE project_options project_warnings OpenMP::OpenMP_CXX)
target_link_libraries(utils PRIVATE project_options project_warnings OpenMP::OpenMP_CXX)

add_executable(Spatial_Distancing main.cpp)
target_link_libraries(Spatial_Distancing PUBLIC
        project_options
        project_warnings
        abm::simulation
        abm::analyser
        abm::basic
        abm::io
        abm::utils
        Boost::filesystem
        OpenMP::OpenMP_CXX)

target_include_directories(Spatial_Distancing PRIVATE ${PROJECT_SOURCE_DIR}/src)


