include(CMakeFindDependencyMacro)

find_dependency(cereal)
find_dependency(Eigen3)
find_dependency(fmt)
find_dependency(GMP)
find_dependency(OpenMP)
find_dependency(ZLIB)

include(${CMAKE_CURRENT_LIST_DIR}/AsmcPrepareDecoding_Runtime.cmake)
