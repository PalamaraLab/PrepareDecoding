// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

// Main PrepareDecoding functions


#ifndef PREPAREDECODING_HPP
#define PREPAREDECODING_HPP

#include "DecodingQuantities.hpp"

namespace asmc {
DecodingQuantities prepareDecoding(std::string_view CSFSFile, std::string_view demographicFile,
                                   std::string_view discretizationFile, int coalescentQuantiles,
                                   int mutationAgeIntervals, std::string_view fileRoot, std::string_view freqFile,
                                   double mutRate, unsigned int samples);

} // namespace asmc

#endif  // PREPAREDECODING_HPP
