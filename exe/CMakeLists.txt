# This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
# See accompanying LICENSE and COPYING for copyright notice and full details.

find_package(cxxopts CONFIG REQUIRED)

add_executable(prepare_decoding_exe prepare_decoding.cpp)
target_link_libraries(prepare_decoding_exe PRIVATE cxxopts::cxxopts)
target_link_libraries(prepare_decoding_exe PRIVATE prepare_decoding_lib project_warnings project_settings)
