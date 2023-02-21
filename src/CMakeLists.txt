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
        DefaultDemographies.cpp
        Utils.cpp
        ArraySpectrum.cpp
        Transition.cpp
        DecodingQuantities.cpp
        PrepareDecoding.cpp
)

set(
        prepare_decoding_hdr
        library.hpp
        ArraySpectrum.hpp
        Csfs.hpp
        CsfsEntry.hpp
        Data.hpp
        DecodingQuantities.hpp
        DefaultArraySpectra.hpp
        DefaultDemographies.hpp
        EigenTypes.hpp
        PrepareDecoding.hpp
        ThinParameterTypes.hpp
        Transition.hpp
        Utils.hpp
)

set(
        prepare_decoding_public_hdr
        ${CMAKE_CURRENT_SOURCE_DIR}/ArraySpectrum.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Csfs.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/CsfsEntry.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Data.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/DecodingQuantities.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/DefaultArraySpectra.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/DefaultDemographies.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/EigenTypes.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/PrepareDecoding.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/ThinParameterTypes.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Transition.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Utils.hpp
)

add_library(prepare_decoding_lib STATIC ${prepare_decoding_src} ${prepare_decoding_hdr})

target_include_directories(
        prepare_decoding_lib INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:>
)
set_target_properties(prepare_decoding_lib PROPERTIES PUBLIC_HEADER "${prepare_decoding_public_hdr}")


target_link_libraries(prepare_decoding_lib PRIVATE Eigen3::Eigen fmt::fmt ZLIB::ZLIB)
target_link_libraries(prepare_decoding_lib PRIVATE project_warnings project_settings)

option(BUILD_WITH_SMCPP_SUPPORT "Enable compilation of smcpp support" ON)
if (BUILD_WITH_SMCPP_SUPPORT)

    find_package(GMP REQUIRED)
    message(STATUS "Found GMP_VERSION:        ${GMP_VERSION}")
    message(STATUS "Found GMP_INCLUDE_DIR:    ${GMP_INCLUDE_DIR}")
    message(STATUS "Found GMP_LIB:            ${GMP_LIB}")
    message(STATUS "Found GMPXX_INCLUDE_DIR:  ${GMPXX_INCLUDE_DIR}")
    message(STATUS "Found GMPXX_LIB:          ${GMPXX_LIB}")

    find_package(MPFR REQUIRED)
    message(STATUS "Found MPFR_FOUND:     ${MPFR_FOUND}")
    message(STATUS "Found MPFR_INCLUDES:  ${MPFR_INCLUDES}")
    message(STATUS "Found MPFR_LIBRARIES: ${MPFR_LIBRARIES}")
    message(STATUS "Found MPFR_VERSION:   ${MPFR_VERSION}")

    if (${MPFR_VERSION} VERSION_LESS 4.0.0)
        add_compile_definitions(MPFR_USE_DEPRECATED_FUNCTIONS)
    endif ()

    find_package(PkgConfig)
    pkg_check_modules(mpfr REQUIRED IMPORTED_TARGET mpfr)

    find_package(cereal CONFIG REQUIRED)
    message(STATUS "Found cereal ${cereal_VERSION}")

    find_package(OpenMP REQUIRED)

    add_library(smcpp_lib 3rd_party/smcpp.cpp 3rd_party/smcpp.hpp 3rd_party/mpreal.h)

    target_include_directories(
            smcpp_lib INTERFACE
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/3rd_party>
            $<INSTALL_INTERFACE:3rd_party>
    )

    target_link_libraries(smcpp_lib PRIVATE Eigen3::Eigen fmt::fmt GMP::gmp cereal OpenMP::OpenMP_CXX PkgConfig::mpfr)
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
