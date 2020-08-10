# This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
# See accompanying LICENSE and COPYING for copyright notice and full details.

find_package(fmt CONFIG REQUIRED)

set(prepare_decoding_src library.cpp)
set(prepare_decoding_hdr library.hpp)

add_library(prepare_decoding_lib STATIC ${prepare_decoding_src} ${prepare_decoding_hdr})
target_include_directories(prepare_decoding_lib INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(prepare_decoding_lib PRIVATE fmt::fmt)
target_link_libraries(prepare_decoding_lib PRIVATE project_warnings project_settings)
