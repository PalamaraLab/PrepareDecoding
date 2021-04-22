# Release Notes

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