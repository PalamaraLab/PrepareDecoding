from asmc.preparedecoding_python_bindings import DecodingQuantities
from asmc.preparedecoding_python_bindings import prepareDecoding
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
            return Discretization(discretization, additional)
        else:
            print("Invalid discretization: expected a path to a file, or a list of the form [[a, b], [c, d], e], where"
                  "each tuple [a, b] is a number and quantity, e.g. [15.0, 2] and e is an optional additional number"
                  "of quantiles to be calculated.")
            sys.exit(1)


def prepare_decoding(
        demography: str,
        discretization: Union[list, str],
        frequencies: str,
        csfs_file: str = "",
        file_root: str = "",
        samples: int = DEFAULT_SAMPLES,
        mutation_rate: float = DEFAULT_MU,

) -> DecodingQuantities:
    """
    Calculate decoding quantities. If a csfs_file is specified, the precalculated CSFS will be used. If no csfs_file is
    specified, CSFS will be calculated.

    :param str demography: the demographic file or code (e.g. 'CEU'),
        NOTE: Effective population size is specified as haploid.
    :param discretization: a discretization file or discretization quantile information,
        for example `[[30.0, 10], [100.0, 8], 25]` specifies a discretization with 
        10 intervals of size 30 generations, followed by 8 intervals of size 100, 
        followed by 25 additional quantiles from the coalescent distribution
    :param str frequencies: the frequencies file, or built-in (e.g. 'UKBB')
    :param str csfs_file: optional file containing precalculated CSFS (default, CSFS will be calculated at runtime)
    :param str file_root: optional file root containing data from which frequencies may be calculated
    :param int samples: number of haploid samples (default: 300)
    :param float mutation_rate: the mutation rate per basepair per generation (default: 1.65e-8)
    :return: a decoding quantities object
    """
    return prepareDecoding(
        demography=Demography(demography),
        discretization=_validate_discretization(discretization),
        frequencies=Frequencies(frequencies, samples),
        csfs_file=csfs_file,
        file_root=file_root,
        samples=samples,
        mut_rate=mutation_rate,
    )
