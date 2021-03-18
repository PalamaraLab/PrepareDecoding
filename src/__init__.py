import sys
import numpy as np

from asmc.preparedecoding_python_bindings import DecodingQuantities
from asmc.preparedecoding_python_bindings import prepareDecoding
from asmc.preparedecoding_python_bindings import CSFS
from asmc.preparedecoding_python_bindings import VectorDouble
from asmc.preparedecoding_python_bindings import VectorEigenMatrix

DEFAULT_MU = 1.65e-8
DEFAULT_SAMPLES = 300


def _make_csfs(
        demographic_file: str,
        discretization_file: str,
        samples: int,
        mu: float = DEFAULT_MU,
) -> CSFS:
    """Make CSFS object using smcpp"""

    try:
        from smcpp import _smcpp
        from smcpp.model import OldStyleModel
    except ImportError:
        error_massage = """This method requires PrepareDecoding be built with optional smcpp dependency.
This (smcpp) is not available on PyPI, so it cannot be installed automatically
when installing PrepareDecoding from PyPI. If you want to generate CSFS from
a demographic file and discretization file, you must also install smcpp:
```
python -m pip install git+https://github.com/popgenmethods/smcpp/@v1.15.3
```
This will require several other dependencies to also be installed. See
https://github.com/PalamaraLab/PrepareDecoding/blob/master/README.md
for details.

You can still create decoding quantities with pre-computed CSFS without
installing smcpp.
"""
        print(error_massage, file=sys.stderr)
        sys.exit(1)

    demo = np.loadtxt(demographic_file)
    array_time = demo[:, 0]
    array_size = demo[:, 1]
    array_disc = np.loadtxt(discretization_file)

    # add dummy last time to get np.diff
    array_time_append = np.append(array_time, array_time[-1] + 100)

    # set n and N0
    n = samples  # number of total haploids, distinguished+undistinguished
    n0 = array_size[0]
    # Population scaled mutation rate
    theta = mu * 2.0 * n0
    om = OldStyleModel(
        array_size / (2.0 * n0),
        array_size / (2.0 * n0),
        np.diff(array_time_append / (2.0 * n0)),
        n0,
    )

    def _csfs(t0, t1):
        res = _smcpp.raw_sfs(om, n - 2, t0 / (2.0 * n0), t1 / (2.0 * n0)) * theta
        res[0, 0] = 1 - np.sum(res)
        return res

    array_disc_original = array_disc
    array_disc = array_disc / (2.0 * n0)

    array_disc = np.append(array_disc, np.inf)
    array_disc_original = np.append(array_disc_original, np.inf)
    froms = array_disc_original[:-1]
    tos = array_disc_original[1:]
    csfses = [_csfs(t0, t1) for t0, t1 in zip(froms, tos)]
    return CSFS.load(
        VectorDouble(array_time),
        VectorDouble(array_size),
        mu,
        samples,
        VectorDouble(froms),
        VectorDouble(tos),
        VectorEigenMatrix(csfses),
    )


def create_from_scratch(
        demographic_file: str,
        discretization_file: str,
        freq_file: str,
        samples: int = DEFAULT_SAMPLES,
        mu: float = DEFAULT_MU,
) -> DecodingQuantities:
    """Run prepareDecoding and construct CSFS automatically using smcpp"""
    return prepareDecoding(
        _make_csfs(demographic_file, discretization_file, samples, mu),
        demographic_file,
        discretization_file,
        freqFile=freq_file,
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
