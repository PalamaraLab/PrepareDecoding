# This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
# See accompanying LICENSE and COPYING for copyright notice and full details.

# Path to the base directory of the repo, no matter where this script is invoked from
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../" >/dev/null 2>&1 && pwd)"

# Read the list of dependencies from file
dependencies=$(cat "${repo_dir}"/scripts/vcpkg_dependencies)

# Location of the vcpkg bootstrap file:
# if it's missing it's probably because the submodule has not been checked out
bootstrap_file="${repo_dir}"/vcpkg/bootstrap-vcpkg.sh

if [ ! -f "${bootstrap_file}" ]; then
    echo "Error: ${bootstrap_file} not found. Did you check out the vcpkg submodule?"
fi

# Location of the vcpkg executable:
# if it's missing it's probably because the bootstrap script has not yet been called
vcpkg_exe=${repo_dir}/vcpkg/vcpkg
if [ ! -f "${vcpkg_exe}" ]; then
    bash "${bootstrap_file}"
fi

# Actually install the dependencies
${vcpkg_exe} install ${dependencies}

# Tell the user how to make use of the toolchain file
toolchain_file=${repo_dir}/vcpkg/scripts/buildsystems/vcpkg.cmake
echo "Use the toolchain file as follows when configuring your CMake build:"
echo "$ cmake -DCMAKE_TOOLCHAIN_FILE=${toolchain_file} /path/to/PrepareDecoding"
