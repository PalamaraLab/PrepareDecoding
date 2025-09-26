[![Unit tests: Ubuntu](https://github.com/PalamaraLab/PrepareDecoding/workflows/Unit%20tests:%20Ubuntu/badge.svg)](https://github.com/PalamaraLab/PrepareDecoding/actions)
[![Unit tests: macOS](https://github.com/PalamaraLab/PrepareDecoding/workflows/Unit%20tests:%20macOS/badge.svg)](https://github.com/PalamaraLab/PrepareDecoding/actions)
[![Regression test](https://github.com/PalamaraLab/PrepareDecoding/workflows/Regression%20test/badge.svg)](https://github.com/PalamaraLab/PrepareDecoding/actions)

[![Static analysis checks](https://github.com/PalamaraLab/PrepareDecoding/workflows/Static%20analysis%20checks/badge.svg)](https://github.com/PalamaraLab/PrepareDecoding/actions)
[![Sanitiser checks](https://github.com/PalamaraLab/PrepareDecoding/workflows/Sanitiser%20checks/badge.svg)](https://github.com/PalamaraLab/PrepareDecoding/actions)
[![codecov](https://codecov.io/gh/PalamaraLab/PrepareDecoding/branch/master/graph/badge.svg)](https://codecov.io/gh/PalamaraLab/PrepareDecoding)

# ASMC Prepare Decoding

Tool to compute decoding quantities.

This is a Python-wrapped C++ project with prebuilt CPython wheels targeting manylinux_2_28 on Linux and macOS 11.0+.

| Platform \ CPython          | ≤3.8 | 3.9 | 3.10 | 3.11 | 3.12 | 3.13 | 3.14 |
|-----------------------------| ---- | --- | ---- | ---- | ---- | ---- | ---- |
| Linux x86_64                | ❌    | ✅   | ✅    | ✅    | ✅    | ✅    | ✅    |
| Linux aarch64               | ❌    | ✅   | ✅    | ✅    | ✅    | ✅    | ✅    |
| macOS Intel (x86_64)        | ❌    | ✅   | ✅    | ✅    | ✅    | ✅    | ✅    |
| macOS Apple Silicon (arm64) | ❌    | ✅   | ✅    | ✅    | ✅    | ✅    | ✅    |


## Quickstart

### Install the Python package from PyPI

Most functionality is available through a Python package which can be installed with:

```bash
pip install asmc-preparedecoding
```

### Example notebook

Examples for using the Python package can be found in the following Jupyter notebook:
- [creating decoding quantities](https://github.com/PalamaraLab/PrepareDecoding/blob/main/notebooks/CreatingDecodingQuantities.ipynb)

Please note that to run the notebook you should first clone the repository and install Jupyter:

```bash
git clone https://github.com/PalamaraLab/PrepareDecoding.git
cd PrepareDecoding

pip install jupyter
jupyter-notebook notebooks/CreatingDecodingQuantities.ipynb
```

### API documentation

A description of the API can be found here:
- [api docs](https://github.com/PalamaraLab/PrepareDecoding/blob/master/docs/api.md)

### File formats

Descriptions of the file formats used can be found here:
- [file formats](https://github.com/PalamaraLab/PrepareDecoding/blob/master/docs/file_formats.md)

## License

This project is currently released under the GNU General Public License Version 3.
