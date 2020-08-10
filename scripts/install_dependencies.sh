# This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
# See accompanying LICENSE and COPYING for copyright notice and full details.

# List of dependencies to install:
#  - catch2     unit testing framework
#  - cxxopts    lightweight command line parser (alternative to boost::program_options
#  - fmt        fast and safe formatting library, partially standardised in C++20
dependencies="catch2 cxxopts fmt"

# Path to the script, no matter where it's invoked from
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"

# Location of the vcpkg bootstrap file:
# if it's missing it's probably because the submodule has not been checked out
bootstrap_file=$(readlink -m ${script_dir}/../vcpkg/bootstrap-vcpkg.sh)
if [ ! -f ${bootstrap_file} ]; then
    echo "Error: ${bootstrap_file} not found. Did you check out the vcpkg submodule?"
fi

# Location of the vcpkg executable:
# if it's missing it's probably because the bootstrap script has not yet been called
vcpkg_exe=$(readlink -m ${script_dir}/../vcpkg/vcpkg)
if [ ! -f ${vcpkg_exe} ]; then
    bash ${bootstrap_file}
fi

# Actually install the dependencies
${vcpkg_exe} install ${dependencies}

# Tell the user how to make use of the toolchain file
toolchain_file=$(readlink -m ${script_dir}/../vcpkg/scripts/buildsystems/vcpkg.cmake)
echo "Use the toolchain file as follows when configuring your CMake build:"
echo "$ cmake -DCMAKE_TOOLCHAIN_FILE=${toolchain_file} /path/to/PrepareDecoding"
