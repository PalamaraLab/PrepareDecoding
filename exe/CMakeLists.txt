# This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
# See accompanying LICENSE and COPYING for copyright notice and full details.

find_package(cxxopts CONFIG REQUIRED)
message(STATUS "Found cxxopts ${cxxopts_VERSION}")

add_executable(prepare_decoding_exe prepare_decoding.cpp)
target_link_libraries(prepare_decoding_exe PRIVATE cxxopts::cxxopts)
target_link_libraries(prepare_decoding_exe PRIVATE prepare_decoding_lib project_warnings project_settings)
if (PYTHON_BINDINGS)
  find_package(pybind11 CONFIG REQUIRED)
  find_package(Eigen3 CONFIG REQUIRED)
  pybind11_add_module(ASMCPrepareDecoding ../src/pybind.cpp)
  target_link_libraries(ASMCPrepareDecoding PRIVATE Eigen3::Eigen prepare_decoding_lib)
endif ()
