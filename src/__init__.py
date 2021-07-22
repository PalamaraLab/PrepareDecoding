from asmc.preparedecoding_python_bindings import DecodingQuantities
from asmc.preparedecoding_python_bindings import prepareDecodingPrecalculatedCsfs
from asmc.preparedecoding_python_bindings import calculateCsfsAndPrepareDecoding
from asmc.preparedecoding_python_bindings import Demography
from asmc.preparedecoding_python_bindings import Discretization
from asmc.preparedecoding_python_bindings import Frequencies
from asmc.preparedecoding_python_bindings import save_demography

from typing import Iterable, Union
import numbers
import sys

DEFAULT_MU = 1.65e-8
DEFAULT_SAMPLES = 300


def _validate_discretization(discretization):
    if isinstance(discretization, str):
        return Discretization(discretization)
    else:
        valid = len(discretization) > 0
        if isinstance(discretization[-1], numbers.Integral):
            additional = discretization.pop()
        else:
            additional = 0
        for x in discretization:
            if not isinstance(x, list) and not len(x) == 2:
                valid = False
                break
            if not isinstance(x[0], numbers.Real) or not isinstance(x[1], numbers.Integral):
                valid = False
        if valid:
            print(f"##### {discretization}, {additional}")
            return Discretization(discretization, additional)
        else:
            print("Invalid discretization: expected a path to a file, or a list of the form [[a, b], [c, d], e], where"
                  "each tuple [a, b] is a number and quantity, e.g. [15.0, 2] and e is an optional additional number"
                  "of quantiles to be calculated.")
            sys.exit(1)


def calculate_csfs_and_prepare_decoding(
        demography: str,
        discretization: Union[list, str],
        freq_file: str,
        samples: int = DEFAULT_SAMPLES,
        mutation_rate: float = DEFAULT_MU,
) -> DecodingQuantities:
    """
    Compute CSFS values and use those to create decoding quantities.

    :param demography: the demographic file or code (e.g. 'CEU')
    :param discretization: the discretization file or discretization quantile information
    :param freq_file: the frequencies file
    :param samples: number of samples (default 300)
    :param mutation_rate: the mutation rate (default 1.65e-8)
    :return: a decoding quantities object
    """
    disc = _validate_discretization(discretization)

    return calculateCsfsAndPrepareDecoding(
        demography=Demography(demography),
        discretization=disc,
        freqFile=Frequencies(freq_file, samples),
        samples=samples,
        mutRate=mutation_rate,
    )


def prepare_decoding_precalculated_csfs(
        csfs_file: str,
        demography: str,
        discretization: Union[list, str],
        freq_file: str,
        samples: int = DEFAULT_SAMPLES,
        mutation_rate: float = DEFAULT_MU,
) -> DecodingQuantities:
    """
    Create decoding quantities from precomputed CSFS values.

    :param csfs_file: file containing the precomputed CSFS values
    :param demography: the demographic file
    :param discretization: the discretization file
    :param freq_file: the frequencies file
    :param samples: number of samples (default 300)
    :param mutation_rate: the mutation rate (default 1.65e-8)
    :return: a decoding quantities object
    """
    return prepareDecodingPrecalculatedCsfs(
        CSFSFile=csfs_file,
        demography=Demography(demography),
        discretization=_validate_discretization(discretization),
        freqFile=Frequencies(freq_file, samples),
        samples=samples,
        mutRate=mutation_rate,
    )
