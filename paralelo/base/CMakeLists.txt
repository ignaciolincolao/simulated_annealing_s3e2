cmake_minimum_required(VERSION 3.25.2)

if (WIN32)
    find_package(CUDAToolkit)
    set(CMAKE_CUDA_ARCHITECTURES 52 60 61 75 86)
endif (WIN32)

if (UNIX)
    set(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)
    set(CMAKE_CUDA_ARCHITECTURES 52 60 61 75 86 89)
endif (UNIX)

project(paralelo LANGUAGES CUDA CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 20)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/../save)

option(ENABLE_SAVE_DATA "Habilita el RecordManager" ON)
if(ENABLE_SAVE_DATA)
    add_compile_options(-DSAVE_DATA=1)
else()
    add_compile_options(-DSAVE_DATA=0)
endif()

include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB_RECURSE all_SRCS
    "${PROJECT_SOURCE_DIR}/src/*.cpp"
    "${PROJECT_SOURCE_DIR}/src/*.c"
    "${PROJECT_SOURCE_DIR}/src/*.cu"

   
)
list(REMOVE_ITEM all_SRCS "${PROJECT_SOURCE_DIR}/src/main.cpp")
list(REMOVE_ITEM all_SRCS "${PROJECT_SOURCE_DIR}/src/test/grid.cpp")

option(BUILD_SHARED_LIB "Build Shared Library" OFF)

if (BUILD_SHARED_LIB)
    add_library(paralelo SHARED src/main.cpp ${all_SRCS})
    install(TARGETS paralelo LIBRARY DESTINATION /usr/local/lib)
    install(DIRECTORY include/
            DESTINATION /usr/local/include/paralelo
            FILES_MATCHING 
            PATTERN "*.h"
            PATTERN "*.hpp"
            PATTERN "*.cuh")
else()
    add_executable(paralelo src/main.cpp ${all_SRCS} alumnos_utm.txt colegios_utm.txt)
    add_executable(grid_executable src/test/grid.cpp ${all_SRCS}  alumnos_utm.txt colegios_utm.txt)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/alumnos_utm.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/colegios_utm.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

    
    set_target_properties(grid_executable PROPERTIES LINKER_LANGUAGE CXX)
    set_property(TARGET grid_executable PROPERTY CUDA_ARCHITECTURES native)
    set_property(TARGET grid_executable PROPERTY CUDA_SEPARABLE_COMPILATION ON)


    target_compile_options(grid_executable PRIVATE $<$<COMPILE_LANGUAGE:CUDA>: 
        -lineinfo --default-stream per-thread 
        --expt-relaxed-constexpr
        -maxrregcount=64
        -Xptxas 
        -dlcm=cg
        -ftz=true 
        -prec-div=false 
        -prec-sqrt=false
    >) 
endif()

set_target_properties(paralelo PROPERTIES LINKER_LANGUAGE CXX)
set_property(TARGET paralelo PROPERTY CUDA_ARCHITECTURES native)
set_property(TARGET paralelo PROPERTY CUDA_SEPARABLE_COMPILATION ON)




target_compile_options(paralelo PRIVATE $<$<COMPILE_LANGUAGE:CUDA>: 
    -lineinfo --default-stream per-thread 
    --expt-relaxed-constexpr
    -maxrregcount=64
    -Xptxas 
    -dlcm=cg
    -ftz=true 
    -prec-div=false 
    -prec-sqrt=false
>) 

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

include(CPack)
