:: This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
:: See accompanying LICENSE and COPYING for copyright notice and full details.

@echo off

:: Path to the base directory of the repo, no matter where this script is invoked from
for %%B in (%~dp0\.) do set repo_dir=%%~dpB

:: Read the list of dependencies from file
set /p dependencies=<%repo_dir%scripts\vcpkg_dependencies

:: Location of the vcpkg bootstrap file:
:: if it's missing it's probably because the submodule has not been checked out
set bootstrap_file=%repo_dir%vcpkg\bootstrap-vcpkg.bat

if not exist %bootstrap_file% (
    echo Error: %bootstrap_file% not found. Did you check out the vcpkg submodule?
    exit /b 1
)

:: Location of the vcpkg executable:
:: if it's missing it's probably because the bootstrap script has not yet been called
set vcpkg_exe=%repo_dir%vcpkg\vcpkg.exe

if not exist %vcpkg_exe% (
    call %bootstrap_file%
)

:: Actually install the dependencies
%vcpkg_exe% install %dependencies%

:: Tell the user how to make use of the toolchain file
set toolchain_file=%repo_dir%vcpkg\scripts\buildsystems\vcpkg.cmake
echo Use the toolchain file as follows when configuring your CMake build:
echo $ cmake -DCMAKE_TOOLCHAIN_FILE=%toolchain_file% path\to\PrepareDecoding
