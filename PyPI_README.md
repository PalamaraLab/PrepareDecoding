[![Unit tests: Windows](https://github.com/PalamaraLab/PrepareDecoding/workflows/Unit%20tests:%20Windows/badge.svg)](https://github.com/PalamaraLab/PrepareDecoding/actions)
[![Unit tests: Ubuntu](https://github.com/PalamaraLab/PrepareDecoding/workflows/Unit%20tests:%20Ubuntu/badge.svg)](https://github.com/PalamaraLab/PrepareDecoding/actions)
[![Unit tests: macOS](https://github.com/PalamaraLab/PrepareDecoding/workflows/Unit%20tests:%20macOS/badge.svg)](https://github.com/PalamaraLab/PrepareDecoding/actions)
[![Regression test](https://github.com/PalamaraLab/PrepareDecoding/workflows/Regression%20test/badge.svg)](https://github.com/PalamaraLab/PrepareDecoding/actions)

[![Static analysis checks](https://github.com/PalamaraLab/PrepareDecoding/workflows/Static%20analysis%20checks/badge.svg)](https://github.com/PalamaraLab/PrepareDecoding/actions)
[![Sanitiser checks](https://github.com/PalamaraLab/PrepareDecoding/workflows/Sanitiser%20checks/badge.svg)](https://github.com/PalamaraLab/PrepareDecoding/actions)
[![codecov](https://codecov.io/gh/PalamaraLab/PrepareDecoding/branch/master/graph/badge.svg)](https://codecov.io/gh/PalamaraLab/PrepareDecoding)
[![BCH compliance](https://bettercodehub.com/edge/badge/PalamaraLab/PrepareDecoding?branch=master)](https://bettercodehub.com/results/PalamaraLab/PrepareDecoding)

# ASMC Prepare Decoding

Tool to compute decoding quantities.

## Quickstart

### Install the Python module from PyPI

Most functionality is available through a Python module which can be installed with:

```bash
pip install asmc-preparedecoding
```

This Python module is currently available on Linux and macOS.
We hope it will be available soon on Windows.

### Example notebook

Examples for using the Python module can be found in the following Jupyter notebook:
- [creating decoding quantities](https://github.com/PalamaraLab/PrepareDecoding/blob/dc870d8a4077498e5c0b35f5a06faa6fdc006422/notebooks/CreatingDecodingQuantities.ipynb)

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
