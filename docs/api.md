# Prepare Decoding API

Import everything from the module:

## Main `prepare_decoding` method

The module `asmc.preparedecoding` contains a method `prepare_decoding` with the following signature:

```python
def prepare_decoding(
        demography: str,
        discretization: Union[list, str],
        frequencies: str,
        csfs_file: str = "",
        samples: int = 300,
        mutation_rate: float = 1.65e-8,
) -> DecodingQuantities
```

### demography

The demography is either a path to a demography file (see [file formats](./file_formats.md)) or a three-letter code representing a built-in demography.

Supported built-in demographies can be one of:
> ACB, ASW, BEB, CDX, CEU, CHB, CHS, CLM, ESN, FIN, GBR, GIH, GWD, IBS, ITU, JPT, KHV, LWK, MSL, MXL, PEL, PJL, PUR, STU, TSI, YRI.

### discretization

The discretization is either a path to a discretization file (see [file formats](./file_formats.md)) or a discretization specification.

A discretization specification allows you to build a discretization with fixed values for recent times with optional additional values calculated at runtime from the coalescent distribution.
This is specified as a list of lists, with an optional final integer, e.g.:

- `[[30.0, 10]]` will create a discretization with 10 intervals of size 30
- `[[30.0, 10], [100.0, 8]]` will create a discretization with 10 intervals of size 30, followed by 8 intervals of size 100
- `[[30.0, 10], [100.0, 8], 25]` will create a discretization with 10 intervals of size 30, followed by 8 intervals of size 100, followed by 25 additional intervals calculated from the coalescent distribution

### frequencies

The frequencies parameter is either a path to a frequencies file (see [file formats](./file_formats.md)) or a built-in code.

Currently, the only built-in available is 'UKBB', which uses built-in frequency information from the UK BioBank, if and only if the number of samples is one of {50, 100, 200, 300}.

### csfs file

An optional path to a CSFS file.
If omitted, the CSFS will be calculated at runtime.

### samples

The number of samples, with a default value of 300.

### mutation rate

The mutation rate, with a default value of 1.65e-8.


## The `DecodingQuantities` object

Calling `prepare_decoding` returns a `DecodingQuantities` object.
This can be explored, if you wish; see the [example notebook](https://github.com/PalamaraLab/PrepareDecoding/blob/dc870d8a4077498e5c0b35f5a06faa6fdc006422/notebooks/CreatingDecodingQuantities.ipynb).

Methods on this object allow you to save information to various files.
For a `DecodingQuantities` object `dq` you can:

- Create a file `dir/output.decodingQuantities.gz`:
    ```python
    dq.save_decoding_quantities('dir/output')
    ```

- Create a file `dir/output.csfs`:
    ```python
    dq.save_csfs('dir/output')
    ```

- Create a file `dir/output.intervalsInfo`:
    ```python
    dq.save_intervals('dir/output')
    ```

- Create a file `dir/output.disc`:
    ```python
    dq.save_discretization('dir/output')
    ```

## Other methods

You can save an in-built demography to file:

- Create file `dir/CEU.demo`:
    ```python
    save_demography('dir/', 'CEU')
    ```
