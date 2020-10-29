import numpy as np
from smcpp import _smcpp, util, estimation_tools
from smcpp.model import OldStyleModel

from ..ASMCPrepareDecoding import (
    DecodingQuantities,
    prepareDecoding,
    CSFSEntry,
    CSFS,
    VectorDouble,
    VectorEigenMatrix,
)

DEFAULT_MU = 1.65e-8
DEFAULT_SAMPLES = 300


def makeCSFS(
    demographicFile: str, discretizationFile: str, samples: int, mu: float = DEFAULT_MU
) -> CSFS:
    "Make CSFS object using smcpp"
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
        samples,
        mu,
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
    "Run prepareDecoding and construct CSFS automatically using smcpp"
    return prepareDecoding(
        makeCSFS(demographicFile, discretizationFile, samples, mu),
        demographicFile,
        discretizationFile,
        freqFile=freqFile,
        samples=samples,
    )
