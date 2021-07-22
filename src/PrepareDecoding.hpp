// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

// Main PrepareDecoding functions


#ifndef PREPAREDECODING_HPP
#define PREPAREDECODING_HPP

#include "DecodingQuantities.hpp"
#include "ThinParameterTypes.hpp"

namespace asmc {
DecodingQuantities prepareDecoding(CSFS& csfs, const Demography& demo, const Discretization& disc,
                                   std::string_view fileRoot, const Frequencies& freq, double mutRate,
                                   unsigned int samples, std::vector<double> discValues);

DecodingQuantities prepareDecodingPrecalculatedCsfs(std::string_view CSFSFile, const Demography& demo,
                                                    const Discretization& disc, std::string_view fileRoot,
                                                    const Frequencies& freq, double mutRate, unsigned int samples);

DecodingQuantities calculateCsfsAndPrepareDecoding(const Demography& demo, const Discretization& disc,
                                                   std::string_view fileRoot, const Frequencies& freq, double mutRate,
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

#endif  // PREPAREDECODING_HPP
