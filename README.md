[![Unit tests: Windows](https://github.com/PalamaraLab/PrepareDecoding/workflows/Unit%20tests:%20Windows/badge.svg)](https://github.com/PalamaraLab/PrepareDecoding/actions)
[![Unit tests: Ubuntu](https://github.com/PalamaraLab/PrepareDecoding/workflows/Unit%20tests:%20Ubuntu/badge.svg)](https://github.com/PalamaraLab/PrepareDecoding/actions)
[![Unit tests: macOS](https://github.com/PalamaraLab/PrepareDecoding/workflows/Unit%20tests:%20macOS/badge.svg)](https://github.com/PalamaraLab/PrepareDecoding/actions)

[![Static analysis checks](https://github.com/PalamaraLab/PrepareDecoding/workflows/Static%20analysis%20checks/badge.svg)](https://github.com/PalamaraLab/PrepareDecoding/actions)
[![Sanitiser checks](https://github.com/PalamaraLab/PrepareDecoding/workflows/Sanitiser%20checks/badge.svg)](https://github.com/PalamaraLab/PrepareDecoding/actions)
[![codecov](https://codecov.io/gh/PalamaraLab/PrepareDecoding/branch/master/graph/badge.svg)](https://codecov.io/gh/PalamaraLab/PrepareDecoding)
[![BCH compliance](https://bettercodehub.com/edge/badge/PalamaraLab/PrepareDecoding?branch=master)](https://bettercodehub.com/results/PalamaraLab/PrepareDecoding)

# Prepare Decoding

Tool to compute decoding quantities.

:warning: this is a work in progress, and is not yet fully featured :warning:

## Quickstart

### Obtaining the source code

Get the source, together with the [vcpkg](https://github.com/microsoft/vcpkg) submodule for installing dependencies:

```bash
git clone --recurse-submodules https://github.com/PalamaraLab/PrepareDecoding.git
cd PrepareDecoding
```

One can either compile the prepare decoding library or use `pip` to install the Python module. If installing via pip, switch to the source folder and run:

    pip install .

The python module is called `ASMCPrepareDecoding` ([see
example](notebooks/PrepareDecoding.ipynb)).

If you wish to install the library or the command line tool, continue to the
next section.

### Installing dependencies

The recommended way to install dependencies is to use the provided scripts, which use the [vcpkg](https://github.com/microsoft/vcpkg) submodule.

On **macOS** and **Linux**, run

```bash
scripts/install_dependencies.sh
```

On **Windows**, run

```bash
scripts\install_dependencies.bat
```

### Configuring and compiling the project

This project uses [CMake](https://cmake.org/).

Create a build directory:

```bash
mkdir build
cd build
```

Configure using the vcpkg toolchain file:

```bash
cmake .. -DCMAKE_TOOLCHAIN_FILE=../vcpkg/scripts/buildsystems/vcpkg.cmake
```

Build the PrepareDecoding library:

```bash
cmake --build . --parallel 4 --target prepare_decoding_lib
```

Build the PrepareDecoding executable:

```bash
cmake --build . --parallel 4 --target prepare_decoding_exe
```

Build the and run the unit tests:

```bash
cmake --build . --parallel 4 --target unit_tests
ctest --output-on-failure
```

## Extra tools

### Coverage

Configure with coverage enabled, using `g++` or `clang++`:

```bash
cmake .. -DCMAKE_TOOLCHAIN_FILE=../vcpkg/scripts/buildsystems/vcpkg.cmake -DCMAKE_BUILD_TYPE=Debug -DENABLE_COVERAGE=ON
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
cmake .. -DCMAKE_TOOLCHAIN_FILE=../vcpkg/scripts/buildsystems/vcpkg.cmake -DCMAKE_BUILD_TYPE=Debug -DENABLE_CLANG_TIDY=ON
cmake --build . --parallel 2 --target prepare_decoding_lib prepare_decoding_exe
```

For cppcheck:

```bash
cmake .. -DCMAKE_TOOLCHAIN_FILE=../vcpkg/scripts/buildsystems/vcpkg.cmake -DCMAKE_BUILD_TYPE=Debug -DENABLE_CPPCHECK=ON
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

## License

This project is currently released under the GNU General Public License Version 3.
