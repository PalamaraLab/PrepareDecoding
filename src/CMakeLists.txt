find_package(fmt CONFIG REQUIRED)

set(ASMC_source library.cpp)
set(ASMC_headers library.hpp)

add_library(prepare_decoding STATIC ${ASMC_source} ${ASMC_headers})
target_include_directories(prepare_decoding INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(prepare_decoding PRIVATE fmt::fmt)
target_link_libraries(prepare_decoding PRIVATE project_warnings project_settings)
