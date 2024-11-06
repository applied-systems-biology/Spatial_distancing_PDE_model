add_library(utils SHARED
        io_util.cpp
        misc_util.cpp
        time_util.cpp)
add_library(abm::utils ALIAS utils)
target_include_directories(utils PRIVATE ${PROJECT_SOURCE_DIR}/src)