# This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
# See accompanying LICENSE and COPYING for copyright notice and full details.

cmake_minimum_required(VERSION 3.16)
message(STATUS "Using CMake version ${CMAKE_VERSION}")
include(cmake/AutodetectVcpkgToolchainFile.cmake)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

project(PrepareDecoding LANGUAGES CXX VERSION 2.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Link this 'library' to set the C++ standard / compile-time options requested
include(cmake/ProjectSettings.cmake)
add_library(project_settings INTERFACE)
target_compile_features(project_settings INTERFACE cxx_std_17)

# Link this 'library' to use the warnings specified in CompilerWarnings.cmake
add_library(project_warnings INTERFACE)
include(cmake/CompilerWarnings.cmake)
set_project_warnings(project_warnings)

# Sanitiser options if supported by compiler
include(cmake/Sanitisers.cmake)
enable_sanitisers(project_settings)

# allow for static analysis options
include(cmake/StaticAnalysers.cmake)

option(BUILD_SHARED_LIBS "Enable compilation of shared libraries" OFF)
option(ENABLE_TESTING "Enable Test Builds" ON)

option(PYTHON_BINDINGS "Whether to build the python bindings" OFF)
if (EXISTS ${CMAKE_SOURCE_DIR}/pybind11/)
    set(PYTHON_BINDINGS ON)
endif()
if (PYTHON_BINDINGS)
    add_subdirectory(pybind11)
endif ()

add_subdirectory(src)
add_subdirectory(exe)

if(ENABLE_TESTING)
    enable_testing()
    add_subdirectory(test)
endif()

option(ENABLE_UNITY "Enable Unity builds of projects" OFF)
if(ENABLE_UNITY)
    # Add for any project you want to apply unity builds for
    set_target_properties(intro PROPERTIES UNITY_BUILD ON)
endif()

add_library(asmc::prepare_decoding ALIAS prepare_decoding_lib)
add_library(asmc::project_warnings ALIAS project_warnings)
add_library(asmc::project_settings ALIAS project_settings)
add_library(asmc::smcpp_lib ALIAS smcpp_lib)
add_executable(asmc::prepare_decoding_exe ALIAS prepare_decoding_exe)
include(GNUInstallDirs)
install(TARGETS prepare_decoding_lib project_warnings project_settings smcpp_lib EXPORT asmc)


install(
        EXPORT asmc
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/AsmcPrepareDecoding
        NAMESPACE AsmcPrepareDecoding::
        FILE AsmcPrepareDecoding_Runtime.cmake
        COMPONENT AsmcPrepareDecoding_Runtime
)

install(
        FILES cmake/AsmcPrepareDecodingConfig.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/AsmcPrepareDecoding
)
