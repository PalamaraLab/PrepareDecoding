find_package(Catch2 CONFIG REQUIRED)

add_library(catch_main STATIC catch_main.cpp)
target_link_libraries(catch_main PUBLIC Catch2::Catch2)
target_link_libraries(catch_main PRIVATE project_settings)

add_executable(unit_tests tests.cpp)
target_link_libraries(unit_tests PRIVATE prepare_decoding project_warnings project_settings catch_main)

# automatically discover tests related to the "tests"
include(Catch)
catch_discover_tests(unit_tests)
