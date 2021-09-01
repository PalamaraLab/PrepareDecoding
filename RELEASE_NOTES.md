# Release Notes

## v2.2.1 (2021-09-01)

Very minor fix to links in documentation.
No change in functionality.

## v2.2 (2021-09-01)

You can now specify discretizations in the following manner:
- as a file: `discretization='/path/to/discretization.disc'` (existing functionality)
- as a number of quantiles, which will be calculated at runtime: `discretization=[100]`
- as a number of pre-specified quantiles plus a number of additional quantiles calculated at runtime: `discretization=[[30.0, 12], [100.0, 15], 39]`
  - this will create 12 discretization points at a spacing of 30.0 (starting from 0.0), followed by 15 at a spacing of 100, followed by 39 additional quantiles

You can now specify built-in frequencies information. Currently, the only supported frequencies are from UKBB. Frequencies can now be specified in the following manner:
- as a file: `frequencies='/path/to/frequencies.frq'` (existing functionality)
- as a string: `frequencies='UKBB'`

### Breaking changes

- The Python API has been simplified, and the strong types mentioned in the v2.1 release are no longer required in Python. 
  Please see the Jupyter Notebook for examples of the current API.
- The strong types remain in the C++ API.
  Please see the file `TestPrepareDecoding.cpp` for examples of the C++ library API.
- There is now a single top-level method `prepare_decoding` (Python) and `prepareDecoding` (C++).
  If the CSFS file parameter is a valid file, CSFS will be loaded from file.
  If the CSFS file parameter is an empty string, CSFS will be calculated at runtime.
  
### Other changes

Various other minor changes have been made.

## v2.1 (2021-05-13)

Default demographies are now bundled with Prepare Decoding.
You can now either supply your own demography file, or choose from the following default demographies:
- ACB, ASW, BEB, CDX, CEU, CHB, CHS, CLM, ESN, FIN, GBR, GIH, GWD, IBS, ITU, JPT, KHV, LWK, MSL, MXL, PEL, PJL, PUR, STU, TSI, YRI

### Breaking changes

- When using the C++ or Python libraries, methods that previous specified a demography file as a string now require an instance of a lightweight strong type `Demography`:
  - In C++, to specify a file: `Demography d("/path/to/demography.demo");`
  - In C++, to specify a default: `Demography d("CEU");`
  - In Python, to specify a file: `d = Demography('/path/to/demography.demo')`
  - In C++, to specify a default: `d = Demography('CEU')`
- A default-constructed Demography will use the default CEU.

### Other changes

- None

## v2.0 (2021-04-22)

Computing CSFS values is now bundled with this project, and no longer relies on the optional `smcpp` dependency.

### Breaking changes

- Python method `create_from_precomputed_csfs` is renamed `prepare_decoding_precalculated_csfs`, but it behaves identically.
- Python method `create_from_scratch` is renamed `calculate_csfs_and_prepare_decoding`, and no longer requires the Python package `smcpp`.

### Other changes

- Some floating point numbers in output files are now written with higher precision.

## v1.1 (2021-03-19)

Minor fixes.

### Breaking changes

- None

### Other changes

- Python wheels now also built for Windows as well as macOS and Linux.
- Corrected syntax in package README.

## v1.0 (2021-03-18)

First public release of ASMC Prepare Decoding, with functionality as described and used in [these notebooks](https://github.com/PalamaraLab/PrepareDecoding/tree/master/notebooks).