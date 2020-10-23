# This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
# See accompanying LICENSE and COPYING for copyright notice and full details.

find_package(Eigen3 CONFIG REQUIRED)
message(STATUS "Found Eigen3 ${Eigen3_VERSION}")

find_package(fmt CONFIG REQUIRED)
message(STATUS "Found {fmt} ${fmt_VERSION}")

find_package(ZLIB REQUIRED)
message(STATUS "Found zlib ${ZLIB_VERSION_STRING}")

set(
        prepare_decoding_src
        library.cpp
        Csfs.cpp
        CsfsEntry.cpp
        Data.cpp
        Utils.cpp
        ArraySpectrum.cpp
        Transition.cpp
        DecodingQuantities.cpp
        PrepareDecoding.cpp
)

set(
        prepare_decoding_hdr
        library.hpp
        Csfs.hpp
        CsfsEntry.hpp
        Data.hpp
        EigenTypes.hpp
        Utils.hpp
        ArraySpectrum.hpp
        Transition.hpp
        DecodingQuantities.hpp
        PrepareDecoding.hpp
)

add_library(prepare_decoding_lib STATIC ${prepare_decoding_src} ${prepare_decoding_hdr})
target_include_directories(prepare_decoding_lib INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(prepare_decoding_lib PRIVATE Eigen3::Eigen fmt::fmt ZLIB::ZLIB)
target_link_libraries(prepare_decoding_lib PRIVATE project_warnings project_settings)
