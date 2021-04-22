from asmc.preparedecoding_python_bindings import DecodingQuantities
from asmc.preparedecoding_python_bindings import prepareDecodingPrecalculatedCsfs
from asmc.preparedecoding_python_bindings import calculateCsfsAndPrepareDecoding

DEFAULT_MU = 1.65e-8
DEFAULT_SAMPLES = 300


def calculate_csfs_and_prepare_decoding(
        demographic_file: str,
        discretization_file: str,
        freq_file: str,
        samples: int = DEFAULT_SAMPLES,
        mutation_rate: float = DEFAULT_MU,
) -> DecodingQuantities:
    """
    Compute CSFS values and use those to create decoding quantities.

    :param demographic_file: the demographic file
    :param discretization_file: the discretization file
    :param freq_file: the frequencies file
    :param samples: number of samples (default 300)
    :param mutation_rate: the mutation rate (default 1.65e-8)
    :return: a decoding quantities object
    """
    return calculateCsfsAndPrepareDecoding(
        demographicFile=demographic_file,
        discretizationFile=discretization_file,
        freqFile=freq_file,
        samples=samples,
        mutRate=mutation_rate,
    )


def prepare_decoding_precalculated_csfs(
        csfs_file: str,
        demographic_file: str,
        discretization_file: str,
        freq_file: str,
        samples: int = DEFAULT_SAMPLES,
        mutation_rate: float = DEFAULT_MU,
) -> DecodingQuantities:
    """
    Create decoding quantities from precomputed CSFS values.

    :param csfs_file: file containing the precomputed CSFS values
    :param demographic_file: the demographic file
    :param discretization_file: the discretization file
    :param freq_file: the frequencies file
    :param samples: number of samples (default 300)
    :param mutation_rate: the mutation rate (default 1.65e-8)
    :return: a decoding quantities object
    """
    return prepareDecodingPrecalculatedCsfs(
        CSFSFile=csfs_file,
        demographicFile=demographic_file,
        discretizationFile=discretization_file,
        freqFile=freq_file,
        samples=samples,
        mutRate=mutation_rate,
    )
