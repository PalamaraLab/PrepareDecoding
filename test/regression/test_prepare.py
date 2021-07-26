import pathlib

from asmc.preparedecoding import *

files_dir = (pathlib.Path(__file__).parents[0]).resolve()

demo_file = str(files_dir / 'input_CEU.demo')
disc_file = str(files_dir / 'input_30-100-2000.disc')
freq_file = str(files_dir / 'input_UKBB.frq')


def from_existing_files():
    """
    Calculate CSFS and prepare decoding using existing files
    """
    dq = prepare_decoding(
        demography=demo_file,
        discretization=disc_file,
        frequencies=freq_file,
        samples=50,
    )

    dq.save_csfs(str(files_dir / 'dq_py_original_files'))
    dq.save_intervals(str(files_dir / 'dq_py_original_files'))
    dq.save_decoding_quantities(str(files_dir / 'dq_py_original_files'))


def builtins_except_frequency():
    """
    Calculate CSFS and prepare decoding using builtins, except frequencies which still comes from file
    """
    dq = prepare_decoding(
        demography='CEU',
        discretization=[[30, 15], [100, 15], 40],
        frequencies=freq_file,
        samples=100,
    )

    dq.save_csfs(str(files_dir / 'dq_py_builtin_except_freq'))
    dq.save_intervals(str(files_dir / 'dq_py_builtin_except_freq'))
    dq.save_decoding_quantities(str(files_dir / 'dq_py_builtin_except_freq'))


def builtins_including_frequency():
    """
    Calculate CSFS and prepare decoding using builtins, including frequencies
    """
    dq = prepare_decoding(
        demography='CEU',
        discretization=[[30, 15], [100, 15], 40],
        frequencies='UKBB',
        samples=100,
    )

    dq.save_csfs(str(files_dir / 'dq_py_builtin_inc_freq'))
    dq.save_intervals(str(files_dir / 'dq_py_builtin_inc_freq'))
    dq.save_decoding_quantities(str(files_dir / 'dq_py_builtin_inc_freq'))


if __name__ == '__main__':
    from_existing_files()
    builtins_except_frequency()
    builtins_including_frequency()
