import sys
import numpy as np

from asmc.preparedecoding_python_bindings import (
    DecodingQuantities,
    prepareDecoding,
    CSFS,
    VectorDouble,
    VectorEigenMatrix,
)
# from .decoding_quantities import DecodingQuantities

DEFAULT_MU = 1.65e-8
DEFAULT_SAMPLES = 300


def makeCSFS(
    demographicFile: str, discretizationFile: str, samples: int, mu: float = DEFAULT_MU
) -> CSFS:
    """Make CSFS object using smcpp"""

    try:
        from smcpp import _smcpp
        from smcpp.model import OldStyleModel
    except ImportError:
        error_massage = """This method requires PrepareDecoding be built with optional smcpp dependency.
This (smcpp) is not available on PyPI, so it cannot be installed automatically
when installing PrepareDecoding from PyPI. If you want to generate CSFS from
a demographic file and discretization file, we recommend installing this module
directly from GitHub: https://github.com/PalamaraLab/PrepareDecoding by running
```
pip install .[smcpp]
```
See the repository README for full installation instructions.

You can still create decoding quantities with pre-computed CSFS.
"""
        print(error_massage, file=sys.stderr)
        sys.exit(1)

    demo = np.loadtxt(demographicFile)
    arrayTime = demo[:, 0]
    arraySize = demo[:, 1]
    arrayDisc = np.loadtxt(discretizationFile)

    # add dummy last time to get np.diff
    arrayTimeAppend = np.append(arrayTime, arrayTime[-1] + 100)

    # set n and N0
    n = samples  # number of total haploids, distinguished+undistinguished
    N0 = arraySize[0]
    # Population scaled mutation rate
    theta = mu * 2.0 * N0
    om = OldStyleModel(
        arraySize / (2.0 * N0),
        arraySize / (2.0 * N0),
        np.diff(arrayTimeAppend / (2.0 * N0)),
        N0,
    )

    def _csfs(t0, t1):
        res = _smcpp.raw_sfs(om, n - 2, t0 / (2.0 * N0), t1 / (2.0 * N0)) * theta
        res[0, 0] = 1 - np.sum(res)
        return res

    arrayDiscOriginal = arrayDisc
    arrayDisc = arrayDisc / (2.0 * N0)

    arrayDisc = np.append(arrayDisc, np.inf)
    arrayDiscOriginal = np.append(arrayDiscOriginal, np.inf)
    froms = arrayDiscOriginal[:-1]
    tos = arrayDiscOriginal[1:]
    csfses = [_csfs(t0, t1) for t0, t1 in zip(froms, tos)]
    return CSFS.load(
        VectorDouble(arrayTime),
        VectorDouble(arraySize),
        mu,
        samples,
        VectorDouble(froms),
        VectorDouble(tos),
        VectorEigenMatrix(csfses),
    )


def run(
        demographicFile: str,
        discretizationFile: str,
        freqFile: str,
        samples: int = DEFAULT_SAMPLES,
        mu: float = DEFAULT_MU,
) -> DecodingQuantities:
    """Run prepareDecoding and construct CSFS automatically using smcpp"""
    return prepareDecoding(
        makeCSFS(demographicFile, discretizationFile, samples, mu),
        demographicFile,
        discretizationFile,
        freqFile=freqFile,
        samples=samples,
    )


def create_from_precomputed_csfs(
        csfs_file: str,
        demographic_file: str,
        discretization_file: str,
        freq_file: str,
        samples: int = DEFAULT_SAMPLES
) -> DecodingQuantities:
    """
    Create decoding quantities from precomputed CSFS values.

    :param csfs_file: file containing the precomputed CSFS values
    :param demographic_file: the demographic file
    :param discretization_file: the discretization file
    :param freq_file: the frequencies file
    :param samples: number of samples (default 300)
    :return: a decoding quantities object
    """
    return prepareDecoding(
        CSFS.loadFromFile(csfs_file),
        demographic_file,
        discretization_file,
        freqFile=freq_file,
        samples=samples,
    )
