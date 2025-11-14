cmake_minimum_required(VERSION 3.0.0)
project(s3e2_sas VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories(
        ${PROJECT_SOURCE_DIR}/include
)
enable_testing()

file(GLOB all_SRCS
    "${PROJECT_SOURCE_DIR}/src/*.cpp"
    "${PROJECT_SOURCE_DIR}/src/*.c"
)
add_executable(s3e2_sas ${all_SRCS} alumnos_utm.txt colegios_utm.txt)
set_target_properties(s3e2_sas PROPERTIES LINKER_LANGUAGE CXX)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/alumnos_utm.txt
    ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/colegios_utm.txt
    ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)