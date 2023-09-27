cmake_minimum_required(VERSION 3.25.2)
if (WIN32)
    find_package(CUDAToolkit)
    set(CMAKE_CUDA_ARCHITECTURES 52 60 61 75 86)

endif (WIN32)
if (UNIX)
    
    set(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)
    set(CMAKE_CUDA_ARCHITECTURES 52 60 61 75)
endif (UNIX)

project(paralelo LANGUAGES CUDA CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/../save)



option(ENABLE_SAVE_DATA "Habilita el almacenamiento de registros" ON)
if(ENABLE_SAVE_DATA)
    add_compile_options(-DSAVE_DATA=1)
else()
    add_compile_options(-DSAVE_DATA=0)
endif()



include_directories(
        ${PROJECT_SOURCE_DIR}/include
)
enable_testing()

file(GLOB_RECURSE all_SRCS
    "${PROJECT_SOURCE_DIR}/src/*.cpp"
    "${PROJECT_SOURCE_DIR}/src/*.c"
    "${PROJECT_SOURCE_DIR}/src/*.cu"
)

# add_executable(paralelo ${all_SRCS} alumnos_utm.txt colegios_utm.txt)

add_library(paralelo SHARED ${all_SRCS})



set_target_properties(paralelo PROPERTIES LINKER_LANGUAGE CXX)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/alumnos_utm.txt
    ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/colegios_utm.txt
    ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

    set_property(TARGET paralelo PROPERTY CUDA_SEPARABLE_COMPILATION ON)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
#set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} $ENV{CUDA_NVCC_FLAGS} -lineinfo --default-stream per-thread -ftz=true -prec-div=false -prec-sqrt=false")
target_compile_options(paralelo PRIVATE $<$<COMPILE_LANGUAGE:CUDA>: -lineinfo --default-stream per-thread 
            --expt-relaxed-constexpr
            -maxrregcount=64
            -Xptxas 
            -dlcm=cg
            -ftz=true 
            -prec-div=false 
            -prec-sqrt=false>) 

# Gprof 
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

include(CPack)
