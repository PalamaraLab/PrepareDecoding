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

option(BUILD_WITH_SMCPP_SUPPORT "Enable compilation of smcpp support" ON)
if (BUILD_WITH_SMCPP_SUPPORT)

    find_package(GMP REQUIRED)
    message(STATUS "Found gmp ${GMP_VERSION}")

    find_package(cereal CONFIG REQUIRED)
    message(STATUS "Found cereal ${cereal_VERSION}")

    find_package(OpenMP REQUIRED)

    add_library(smcpp_lib 3rd_party/smcpp.cpp 3rd_party/smcpp.hpp 3rd_party/mpreal.h)
    target_include_directories(smcpp_lib INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party)
    target_link_libraries(smcpp_lib PRIVATE Eigen3::Eigen fmt::fmt gmp gmpxx cereal OpenMP::OpenMP_CXX)
    set_target_properties(smcpp_lib PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
    set_target_properties(smcpp_lib PROPERTIES CXX_CLANG_TIDY "")
    set_target_properties(smcpp_lib PROPERTIES CXX_CPPCHECK "")
    target_link_libraries(prepare_decoding_lib PRIVATE smcpp_lib)
endif ()


if (PYTHON_BINDINGS)
    set_target_properties(prepare_decoding_lib PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
    pybind11_add_module(preparedecoding_python_bindings pybind.cpp)
    target_link_libraries(preparedecoding_python_bindings PRIVATE Eigen3::Eigen fmt::fmt prepare_decoding_lib project_settings)
endif ()
