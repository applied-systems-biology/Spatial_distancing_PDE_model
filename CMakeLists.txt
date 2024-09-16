cmake_minimum_required(VERSION 3.13)

cmake_policy(SET CMP0079 NEW)

project(Spatial_Distancing VERSION 1.0.0)
set(DEBUG_LEVEL "LOG_SYSTEM")

find_package(Boost 1.56 REQUIRED COMPONENTS filesystem)
find_package(OpenMP 4.5 REQUIRED)

add_library(project_options INTERFACE)
target_compile_features(project_options INTERFACE cxx_std_17)
add_compile_definitions(project_options INTERFACE ${DEBUG_LEVEL}=1)

add_library(project_warnings INTERFACE)
include(cmake/CompilerWarnings.cmake)
#set_project_warnings(project_warnings)

# checkout Sanitizers.cmake for enabling sanitizers
include(cmake/Sanitizers.cmake)
enable_sanitizers(project_options)

# checkout Doxygen.cmake for enabling doxygen
include(cmake/Doxygen.cmake)
enable_doxygen()

# allow for static analysis options
include(cmake/StaticAnalyzers.cmake)

option(ENABLE_PCH "Enable Precompiled Headers" OFF)
if(ENABLE_PCH)
    target_precompile_headers(
            project_options
            INTERFACE
            <memory>
            <vector>
            <string>
            <map>
            <utility>)
endif()

add_subdirectory(src)
enable_testing()
add_subdirectory(test)