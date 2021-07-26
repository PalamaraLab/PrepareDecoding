// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

// Main PrepareDecoding functions

#ifndef PREPAREDECODING_HPP
#define PREPAREDECODING_HPP

#include "DecodingQuantities.hpp"
#include "ThinParameterTypes.hpp"

namespace asmc {
DecodingQuantities calculateDecodingQuantities(CSFS& csfs, const Demography& demo, const Discretization& disc,
                                               std::string_view fileRoot, const Frequencies& freq, double mutRate,
                                               unsigned int samples, std::vector<double> discValues);

DecodingQuantities calculateCsfsAndDecodingQuantities(const Demography& demo, const Discretization& disc,
                                                      std::string_view fileRoot, const Frequencies& freq,
                                                      double mutRate, unsigned int samples);

/**
 * Prepare decoding (calculate the decoding quantities).
 *
 * If a valid CSFS File is specified, the pre-calculated CSFS will be used.  Leaving the CSFSFile parameter as an empty
 * string will cause the CSFS to be calculated.
 *
 * @param demo a demography object
 * @param disc a discretization object
 * @param freq a frequencies object
 * @param CSFSFile a path to a CSFS file, or an empty string_view
 * @param fileRoot the data file root, or an empty string
 * @param mutRate the mutation rate
 * @param samples the number of samples
 * @return the decoding quantities calculated from the given parameters
 */
DecodingQuantities prepareDecoding(const Demography& demo, const Discretization& disc, const Frequencies& freq,
                                   std::string_view CSFSFile, std::string_view fileRoot, double mutRate,
                                   unsigned int samples);

/**
 * Read demographic info from file, or use default EU demographic information
 * @param demographicFile path to the demographic file, or an empty stringview
 * @return vector of times, appended with inf, and vector of sizes, appended with a copy of the final element
 */
std::tuple<std::vector<double>, std::vector<double>> getDemographicInfo(const Demography& demo);

/**
 * Read discretization info from file, or use coalescent quantiles or mutation age intervals
 * @param disc
 * @param times
 * @param sizes
 * @return vector of discretizations, appended with inf
 */
std::vector<double> getDiscretizationInfo(const Discretization& disc, const std::vector<double>& times,
                                          const std::vector<double>& sizes);

} // namespace asmc

#endif // PREPAREDECODING_HPP
