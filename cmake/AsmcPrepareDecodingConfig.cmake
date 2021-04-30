include(CMakeFindDependencyMacro)

find_dependency(cereal)
find_dependency(cxxopts)
find_dependency(eigen3)
find_dependency(fmt)
find_dependency(zlib)

include(${CMAKE_CURRENT_LIST_DIR}/AsmcPrepareDecoding_Runtime.cmake)
