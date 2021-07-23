import pathlib

from asmc.preparedecoding import *


# Create a directory alongside this script for files to be written to
files_dir = (pathlib.Path(__file__).parents[0] / 'dq_for_asmc').resolve()
files_dir.mkdir(parents=False, exist_ok=True)


# First generate with a coarse discretization
disc_coarse = [[30.0, 14], [100.0, 15], 39]
dq_coarse = calculate_csfs_and_prepare_decoding(
    demography='CEU',
    discretization=disc_coarse,
    frequencies='UKBB',
    samples=300,
)

dq_coarse.save_csfs(str(files_dir / '30-100-2000_CEU'))
dq_coarse.save_decoding_quantities(str(files_dir / '30-100-2000_CEU'))
dq_coarse.save_intervals(str(files_dir / '30-100-2000_CEU'))
dq_coarse.save_discretization(str(files_dir / '30-100-2000_CEU'))


# Next, generate with a finer discretization
disc_fine = [[10.0, 40], [20.0, 79], 39]
dq_fine = calculate_csfs_and_prepare_decoding(
    demography='CEU',
    discretization=disc_fine,
    frequencies='UKBB',
    samples=300,
)

dq_fine.save_csfs(str(files_dir / '10-20-2000_CEU'))
dq_fine.save_decoding_quantities(str(files_dir / '10-20-2000_CEU'))
dq_fine.save_intervals(str(files_dir / '10-20-2000_CEU'))
dq_fine.save_discretization(str(files_dir / '10-20-2000_CEU'))


# Finally, save the demography
save_demography(str(files_dir), 'CEU')
