[![Unit tests: Windows](https://github.com/PalamaraLab/PrepareDecoding/workflows/Unit%20tests:%20Windows/badge.svg)](https://github.com/PalamaraLab/PrepareDecoding/actions)
[![Unit tests: Ubuntu](https://github.com/PalamaraLab/PrepareDecoding/workflows/Unit%20tests:%20Ubuntu/badge.svg)](https://github.com/PalamaraLab/PrepareDecoding/actions)
[![Unit tests: macOS](https://github.com/PalamaraLab/PrepareDecoding/workflows/Unit%20tests:%20macOS/badge.svg)](https://github.com/PalamaraLab/PrepareDecoding/actions)
[![Regression test](https://github.com/PalamaraLab/PrepareDecoding/workflows/Regression%20test/badge.svg)](https://github.com/PalamaraLab/PrepareDecoding/actions)

[![Static analysis checks](https://github.com/PalamaraLab/PrepareDecoding/workflows/Static%20analysis%20checks/badge.svg)](https://github.com/PalamaraLab/PrepareDecoding/actions)
[![Sanitiser checks](https://github.com/PalamaraLab/PrepareDecoding/workflows/Sanitiser%20checks/badge.svg)](https://github.com/PalamaraLab/PrepareDecoding/actions)
[![codecov](https://codecov.io/gh/PalamaraLab/PrepareDecoding/branch/master/graph/badge.svg)](https://codecov.io/gh/PalamaraLab/PrepareDecoding)
[![BCH compliance](https://bettercodehub.com/edge/badge/PalamaraLab/PrepareDecoding?branch=master)](https://bettercodehub.com/results/PalamaraLab/PrepareDecoding)

# Prepare Decoding

Tool to compute decoding quantities.

## Quickstart

### Install the Python module from PyPI

Most functionality is available through a Python module which can be installed with:

```bash
pip install asmc-preparedecoding
```

This Python module is currently available on Linux and macOS.
We hope it will be available soon on Windows.

Examples for using the Python module can be found in the following Jupyter notebook:
- [creating decoding quantities](notebooks/CreatingDecodingQuantities.ipynb)

Please note that you must install Jupyter in order to view the notebook, and then open it:

```bash
pip install jupyter
jupyter-notebook notebooks/CreatingDecodingQuantities.ipynb
```

### Compiling the C++ library and executable

Get the source, together with the [vcpkg](https://github.com/microsoft/vcpkg) and [pybind11](https://github.com/pybind/pybind11) submodules:

```bash
git clone --recurse-submodules https://github.com/PalamaraLab/PrepareDecoding.git
cd PrepareDecoding
```

The recommended way to install dependencies is via the [vcpkg](https://github.com/microsoft/vcpkg) submodule.
If you have checked out this submodule, most dependencies will be automatically installed when you run the CMake configuration step (below).

Several additional dependencies should be obtained from your package manager:

- Ubuntu/Debian:
    ```bash
    sudo apt install libgmp-dev libmpfr-dev
    ```

- CentOS/Fedora:
    ```bash
    sudo yum install gmp-devel mpfr-devel
    ```

- macOS:
    ```bash
    brew install gmp mpfr libomp
    ```

### Configuring and compiling the project

This project uses [CMake](https://cmake.org/).

Create a build directory:

```bash
mkdir build
cd build
```

Configure, and build the PrepareDecoding library and executable:

```bash
cmake ..
cmake --build . --parallel 4 --target prepare_decoding_lib
cmake --build . --parallel 4 --target prepare_decoding_exe
```

You can optionally build the and run the unit tests:

```bash
cmake --build . --parallel 4 --target unit_tests
ctest --output-on-failure
```

## Extra tools for C++ developers

### Coverage

Configure with coverage enabled, using `g++` or `clang++`:

```bash
cmake .. -DCMAKE_BUILD_TYPE=Debug -DENABLE_COVERAGE=ON
cmake --build . --parallel 2 --target unit_tests
ctest -j2 --output-on-failure
lcov --directory . --capture --output-file coverage.info
lcov --remove coverage.info '/usr/*' '*/vcpkg/*' '*/test/*' --output-file coverage.info
lcov --list coverage.info
```

### Static analysis

This project is configured to work with [clang tidy](https://clang.llvm.org/extra/clang-tidy/) and [cppcheck](http://cppcheck.sourceforge.net/).
Enable the relevant option and compile with `clang++`.

For clang tidy:

```bash
cmake .. -DCMAKE_BUILD_TYPE=Debug -DENABLE_CLANG_TIDY=ON
cmake --build . --parallel 2 --target prepare_decoding_lib prepare_decoding_exe
```

For cppcheck:

```bash
cmake .. -DCMAKE_BUILD_TYPE=Debug -DENABLE_CPPCHECK=ON
cmake --build . --parallel 2 --target prepare_decoding_lib prepare_decoding_exe
```

### Sanitisers

This project is configured to work with various LLVM sanitisers.
Enable the relevant option, compile with `clang++`, and run the unit tests.

```bash
cmake .. -DCMAKE_TOOLCHAIN_FILE=../vcpkg/scripts/buildsystems/vcpkg.cmake -DCMAKE_BUILD_TYPE=Debug -DENABLE_SANITISER_{{{SANITISER}}}=ON
cmake --build . --parallel 2 --target unit_tests
./test/unit_tests
```

where `{{{SANITISER}}}` is one of:

- [ADDRESS](https://clang.llvm.org/docs/AddressSanitizer.html)
- [LEAK](https://clang.llvm.org/docs/LeakSanitizer.html)
- [UNDEFINED_BEHAVIOUR](https://clang.llvm.org/docs/UndefinedBehaviorSanitizer.html)
- [THREAD](https://clang.llvm.org/docs/ThreadSanitizer.html)


## For developers: making a release

- Bump the version number in [setup.py](setup.py), [CMakeLists.txt](CMakeLists.txt), and [vcpkg.json](vcpkg.json)
- Update [RELEASE_NOTES.md](RELEASE_NOTES.md)
- Push changes and check that all [GitHub workflows](https://github.com/PalamaraLab/PrepareDecoding/actions) pass
- Tag the commit in Git using syntax `vX.Y`
- Make a release on GitHub, which should trigger a new build that will upload Python wheels to PyPI


## License

This project is currently released under the GNU General Public License Version 3.
