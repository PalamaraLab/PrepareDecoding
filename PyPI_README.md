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

This Python module is available on Linux, macOS and Windows.

Some additional functionality, for creating CSFS, requires the additional dependency [smcpp](https://github.com/popgenmethods/smcpp/), which is not available via PyPI.
If you require this functionality, you should additionally [follow these instructions](#installing-smcpp) to install smcpp.

Examples for using the Python module can be found in the notebooks directory:
- [creating decoding quantities from precomputed CSFS](https://github.com/PalamaraLab/PrepareDecoding/blob/master/notebooks/CreateDecodingQuantitiesFromPrecomputedCSFS.ipynb)
- [creating decoding quantities from scratch](https://github.com/PalamaraLab/PrepareDecoding/blob/master/notebooks/CreateDecodingQuantitiesFromScratch.ipynb) (requires smcpp)

### Installing smcpp

The optional smcpp dependency is not available on PyPI, and itself requires a few additional dependencies.

On **Linux**, run

```bash
sudo apt install libgmp-dev libmpfr-dev libgsl0-dev
```

on **macOS**, run

```bash
brew install mpfr gmp gsl
```

Then, we recommend starting from a clean virtual environment. 
Switch to the source directory and run:

```bash
python3 -m venv venv
source venv/bin/activate
python -m pip install --upgrade pip setuptools wheel
python -m pip install asmc-preparedecoding
python -m pip install git+https://github.com/popgenmethods/smcpp/@v1.15.3
```

## License

This project is currently released under the GNU General Public License Version 3.
