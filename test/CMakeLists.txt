# This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
# See accompanying LICENSE and COPYING for copyright notice and full details.

Include(FetchContent)
FetchContent_Declare(
        Catch2
        GIT_SHALLOW TRUE
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG 2ab20a0e008845e02bd06248e61ca6e5ad1aba33 # version 3.3.1
)
FetchContent_MakeAvailable(Catch2)
message(STATUS "Obtained Catch2 ${Catch2_VERSION}")

add_definitions(-DPREPARE_DECODING_TEST_DIR=\"${CMAKE_CURRENT_LIST_DIR}\")

set(
        test_src
        tests.cpp
        TestCsfs.cpp
        TestCsfsEntry.cpp
        TestData.cpp
        TestDefaultDemographies.cpp
        TestPrepareDecoding.cpp
        TestThinParameterTypes.cpp
        TestUtils.cpp
)

add_executable(unit_tests ${test_src})
target_link_libraries(unit_tests PRIVATE prepare_decoding_lib project_warnings project_settings Catch2::Catch2WithMain)

# automatically discover tests
list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras)
include(CTest)
include(Catch)
catch_discover_tests(unit_tests)
