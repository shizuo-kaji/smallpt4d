cmake_minimum_required(VERSION 3.16)

project(smallpt4d LANGUAGES CXX)

option(SMALLPT4D_ENABLE_OPENMP "Enable OpenMP parallel rendering." ON)
option(SMALLPT4D_ENABLE_IPO "Enable interprocedural optimization (LTO) in Release." ON)
option(SMALLPT4D_ENABLE_OPENCL "Enable OpenCL GPU backend when available." ON)

add_executable(smallpt4d
    main.cpp
    opencl_renderer.cpp
    scene_definition.cpp
    stdafx.cpp
)

target_compile_features(smallpt4d PRIVATE cxx_std_17)
set_target_properties(smallpt4d PROPERTIES CXX_EXTENSIONS OFF)

if(MSVC)
    target_compile_options(smallpt4d PRIVATE /W4)
else()
    target_compile_options(smallpt4d PRIVATE -Wall)
endif()

# Prefer an installed Eigen3 package; fall back to common header paths.
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(Eigen3_FOUND)
    target_link_libraries(smallpt4d PRIVATE Eigen3::Eigen)
else()
    find_path(EIGEN3_INCLUDE_DIR Eigen/Dense
        HINTS
            ENV EIGEN3_INCLUDE_DIR
            /opt/homebrew/include/eigen3
            /usr/local/include/eigen3
            /usr/include/eigen3
    )

    if(NOT EIGEN3_INCLUDE_DIR)
        message(FATAL_ERROR
            "Could not find Eigen3 headers. Install Eigen3 or set EIGEN3_INCLUDE_DIR."
        )
    endif()

    target_include_directories(smallpt4d PRIVATE "${EIGEN3_INCLUDE_DIR}")
endif()

find_package(Threads REQUIRED)
target_link_libraries(smallpt4d PRIVATE Threads::Threads)

if(SMALLPT4D_ENABLE_IPO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT SMALLPT4D_IPO_SUPPORTED OUTPUT SMALLPT4D_IPO_ERROR)
    if(SMALLPT4D_IPO_SUPPORTED)
        set_target_properties(smallpt4d PROPERTIES INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
    else()
        message(STATUS "IPO/LTO not enabled: ${SMALLPT4D_IPO_ERROR}")
    endif()
endif()

if(SMALLPT4D_ENABLE_OPENMP)
    find_package(OpenMP QUIET)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(smallpt4d PRIVATE OpenMP::OpenMP_CXX)
    else()
        message(WARNING "OpenMP not found; building without OpenMP support.")
    endif()
endif()

if(SMALLPT4D_ENABLE_OPENCL)
    find_package(OpenCL QUIET)
    if(OpenCL_FOUND)
        target_compile_definitions(smallpt4d PRIVATE SMALLPT4D_WITH_OPENCL=1)
        if(TARGET OpenCL::OpenCL)
            target_link_libraries(smallpt4d PRIVATE OpenCL::OpenCL)
        else()
            target_include_directories(smallpt4d PRIVATE ${OpenCL_INCLUDE_DIRS})
            target_link_libraries(smallpt4d PRIVATE ${OpenCL_LIBRARIES})
        endif()
    else()
        message(STATUS "OpenCL not found; GPU backend disabled.")
    endif()
endif()
