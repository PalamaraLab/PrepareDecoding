import os
import re
import sys
import platform
import subprocess

from setuptools import setup, Extension, find_namespace_packages
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(["cmake", "--version"])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: "
                + ", ".join(e.name for e in self.extensions)
            )

        if platform.system() == "Windows":
            cmake_version = LooseVersion(
                re.search(r"version\s*([\d.]+)", out.decode()).group(1)
            )
            if cmake_version < "3.1.0":
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        subprocess.check_call(
            ["bash", os.path.join(ext.sourcedir, "scripts", "install_dependencies.sh")]
        )
        cmake_args = [
            "-DCMAKE_TOOLCHAIN_FILE="
            + os.path.join(ext.sourcedir, "vcpkg/scripts/buildsystems/vcpkg.cmake"),
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + extdir,
            "-DPYTHON_BINDINGS=ON",
            "-DPYTHON_EXECUTABLE=" + sys.executable,
        ]

        cfg = "Debug" if self.debug else "Release"
        warnings_as_errors = "ON" if self.debug else "OFF"

        # build_args = ['--config', cfg]
        build_args = []

        if platform.system() == "Windows":
            cmake_args += [
                "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}".format(cfg.upper(), extdir)
            ]
            if sys.maxsize > 2 ** 32:
                cmake_args += ["-A", "x64"]
            build_args += ["--", "/m"]
        else:
            cmake_args += ["-DCMAKE_BUILD_TYPE=" + cfg]
            cmake_args += ["-DWARNINGS_AS_ERRORS=" + warnings_as_errors]
            build_args += ["--", "-j2"]

        env = os.environ.copy()
        env["CXXFLAGS"] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get("CXXFLAGS", ""), self.distribution.get_version()
        )
        build = os.path.join(self.build_temp, "build")
        if not os.path.exists(build):
            os.makedirs(build)

        print("Using CMake arg " + x for x in cmake_args)

        subprocess.check_call(
            ["cmake", ext.sourcedir] + cmake_args,
            cwd=build,
            env=env,
        )
        subprocess.check_call(
            ["cmake", "--build", "."] + build_args,
            cwd=build,
        )


with open('PyPI_README.md', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="asmc-preparedecoding",
    version="1.0",
    author="PalamaraLab (https://palamaralab.github.io/)",
    install_requires=["numpy"],
    description="Prepare decoding quantities for ASMC",
    packages=find_namespace_packages(include=['asmc.*']),
    long_description=long_description,
    long_description_content_type='text/markdown',
    ext_modules=[CMakeExtension("asmc/preparedecoding")],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
)
