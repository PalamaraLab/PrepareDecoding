// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

// Main PrepareDecoding functions

#include "PrepareDecoding.hpp"
#include "Csfs.hpp"
#include "Data.hpp"
#include "DecodingQuantities.hpp"
#include "Transition.hpp"
#include "Utils.hpp"
#include <cstdlib>
#include <filesystem>
#include <fmt/core.h>

#include "smcpp.hpp"

namespace fs = std::filesystem;

namespace asmc {

class PrepareDecoding {

};


DecodingQuantities prepareDecoding(CSFS& csfs, std::string_view demographicFile, std::string_view discretizationFile,
                                   int coalescentQuantiles, int mutationAgeIntervals, std::string_view fileRoot,
                                   std::string_view freqFile, double mutRate, unsigned int samples) {

  auto [times, sizes] = getDemographicInfo(demographicFile);
  auto discs = getDiscretizationInfo(discretizationFile, coalescentQuantiles, mutationAgeIntervals, times, sizes);



  Data data;
  if (!fileRoot.empty()) {
    fmt::print("Files will be read from: {}*\n", fileRoot);
  }
  if (!freqFile.empty()) {
    if (!fs::exists(freqFile)) {
      throw std::runtime_error(fmt::format("Could not open {}\n", freqFile));
    }
    fmt::print("Will load minor allele frequencies from {} ...\n", freqFile);
    data.addFreq(freqFile);
  } else {
    if (fileRoot.empty()) {
      throw std::runtime_error("Either one of --freqFile or --fileRoot has to be specified\n");
    }
    fmt::print("Will compute allele frequencies in haps files with root {}.", fileRoot);
    data = Data(fileRoot);
  }
  fmt::print("Will use mutation rate mu = {}.\n", mutRate);
  samples = std::min(samples, data.getHaploidSampleSize());
  fmt::print("Number of samples in CSFS calculations: {}.\n", samples);

  // Transition
  Transition transition(times, sizes, discs, TransitionType::CSC);
  if (csfs.verify(times, sizes, mutRate, samples, discs)) {
    fmt::print("Verified " + std::to_string(csfs.getCSFS().size()) + " CSFS entries.\n");
  } else {
    throw std::runtime_error("CSFS could not be verified. The Python version can be used"
                             "to fetch CSFS and prepare decoding quantities.");
  }

  csfs.fixAscertainment(data, samples, transition);

  // Build decoding quantities
  fmt::print("\nBuilding decoding quantities...\n");
  return DecodingQuantities(csfs, transition, mutRate);
}

DecodingQuantities calculateCsfsAndPrepareDecoding(std::string_view demographicFile,
                                                   std::string_view discretizationFile, const int coalescentQuantiles,
                                                   const int mutationAgeIntervals, std::string_view fileRoot,
                                                   std::string_view freqFile, const double mutRate,
                                                   const unsigned int samples) {

  // Get the array times and sizes, and remove the additional element added to the end of each array
  auto [arrayTime, arraySize] = getDemographicInfo(demographicFile);
  arrayTime.pop_back();
  arraySize.pop_back();

  auto arrayDisc =
      getDiscretizationInfo(discretizationFile, coalescentQuantiles, mutationAgeIntervals, arrayTime, arraySize);

  const double N0 = arraySize.front();
  const double theta = mutRate * 2. * N0;

  // Append a value to the end of arrayTime so the finite difference is the same size as the original array
  std::vector arrayTimeAppend = arrayTime;
  arrayTimeAppend.push_back(arrayTimeAppend.back() + 100.0);
  assert(arrayTime.size() + 1ul == arrayTimeAppend.size());

  std::vector<double> aVec;
  std::vector<double> sVec;
  for (auto i = 0ul; i < arraySize.size(); ++i) {
    aVec.push_back(arraySize[i] / (2.0 * N0));
    sVec.push_back((arrayTimeAppend[i + 1] - arrayTimeAppend[i]) / (2.0 * N0));
  }

  std::vector<double> froms;
  std::vector<double> tos;
  std::vector<Matrix<double>> csfses;

  // Must initialise smcpp cache once
  smcpp_init_cache();
  for (auto i = 0ul; i < arrayDisc.size() - 1; ++i) {

    auto t0 = arrayDisc[i];
    auto t1 = arrayDisc[i + 1];

    froms.push_back(t0);
    tos.push_back(t1);

    csfses.emplace_back(raw_sfs(aVec, sVec, static_cast<int>(samples - 2u), t0 / (2. * N0), t1 / (2. * N0)) * theta);
    csfses.back()(0, 0) = 1.0 - csfses.back().sum();
  }

  auto csfs = CSFS::load(arrayTime, arraySize, mutRate, samples, froms, tos, csfses);

  return prepareDecoding(csfs, demographicFile, discretizationFile, coalescentQuantiles, mutationAgeIntervals, fileRoot,
                         freqFile, mutRate, samples);
}

DecodingQuantities prepareDecodingPrecalculatedCsfs(std::string_view CSFSFile, std::string_view demographicFile,
                                           std::string_view discretizationFile, int coalescentQuantiles,
                                           int mutationAgeIntervals, std::string_view fileRoot,
                                           std::string_view freqFile, double mutRate, unsigned int samples) {
  if (!CSFSFile.empty() & fs::exists(CSFSFile)) {
    fmt::print("Will load precomputed CSFS from {} ...\n", CSFSFile);
    auto csfs = CSFS::loadFromFile(CSFSFile);
    return prepareDecoding(csfs, demographicFile, discretizationFile, coalescentQuantiles, mutationAgeIntervals,
                           fileRoot, freqFile, mutRate, samples);
  } else {
    throw std::runtime_error("Valid CSFS file needs to be specified\n");
  }
}

std::tuple<std::vector<double>, std::vector<double>> getDemographicInfo(std::string_view demographicFile) {
  std::vector<double> times;
  std::vector<double> sizes;

  if (!demographicFile.empty() & fs::exists(demographicFile)) {
    fmt::print("Will read demographic model from {} ...\n", demographicFile);
    auto ts = readDemographic(demographicFile);
    times = ts.first;
    sizes = ts.second;
  } else {
    times = std::vector<double>(Transition::EUtime.begin(), Transition::EUtime.end());
    sizes = std::vector<double>(Transition::EUsize.begin(), Transition::EUsize.end());
    fmt::print("Did not input a demographic model, using default EU model.\n");
  }
  times.emplace_back(std::numeric_limits<double>::infinity());
  sizes.emplace_back(sizes.back()); // duplicate last element
  return std::make_tuple(times, sizes);
}

std::vector<double> getDiscretizationInfo(std::string_view discretizationFile, const int coalescentQuantiles,
                                          const int mutationAgeIntervals, const std::vector<double>& times,
                                          const std::vector<double>& sizes) {
  std::vector<double> discs;

  if (!discretizationFile.empty() & fs::exists(discretizationFile)) {
    fmt::print("Will read discretization intervals from {} ...\n", discretizationFile);
    discs = readDiscretization(discretizationFile);
  } else if (coalescentQuantiles > 0) {
    discs = Transition::getTimeExponentialQuantiles(coalescentQuantiles, times, sizes);
    fmt::print("Using {} discretization intervals from coalescent distribution.\n", coalescentQuantiles);
  } else if (mutationAgeIntervals > 0) {
    discs = Transition::getTimeErlangQuantiles(mutationAgeIntervals, times, sizes);
    fmt::print("Using {} discretization intervals from mutation age intervals.\n", mutationAgeIntervals);
  }
  if (discs.empty()) {
    throw std::runtime_error(
        "Specify a valid option from --discretization, --coalescentQuantiles, --mutationQuantiles\n");
  }
  discs.emplace_back(std::numeric_limits<double>::infinity());

  return discs;
}

} // namespace asmc
