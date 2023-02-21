// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

// Main PrepareDecoding functions

#include "PrepareDecoding.hpp"
#include "Csfs.hpp"
#include "Data.hpp"
#include "DecodingQuantities.hpp"
#include "DefaultDemographies.hpp"
#include "Transition.hpp"
#include "Utils.hpp"
#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fmt/core.h>
#include <utility>

#include "smcpp.hpp"

namespace fs = std::filesystem;

namespace asmc {

DecodingQuantities calculateDecodingQuantities(CSFS& csfs, const Demography& demo, const Discretization& disc,
                                               std::string_view fileRoot, const Frequencies& freq, double mutRate,
                                               unsigned int samples, std::vector<double> discValues = {}) {

  auto [times, sizes] = getDemographicInfo(demo);

  if (discValues.empty()) {
    discValues = getDiscretizationInfo(disc, times, sizes);
  }

  Data data;
  if (!fileRoot.empty()) {
    fmt::print("Files will be read from: {}*\n", fileRoot);
  }

  if (freq.isBuiltIn()) {
    assert(freq.getNumSamples() == samples);
    data.addFreq(freq);
    fmt::print("Using built-in frequency information from {} ...\n", freq.getFreqIdentifier());
  } else if (freq.isFile()) {
    fmt::print("Will use minor allele frequencies from {} ...\n", freq.getFreqIdentifier());
    data.addFreq(freq);
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
  Transition transition(times, sizes, discValues, TransitionType::CSC);
  if (csfs.verify(times, sizes, mutRate, samples, discValues)) {
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

DecodingQuantities prepareDecoding(const Demography& demo, const Discretization& disc, const Frequencies& freq,
                                   std::string_view CSFSFile, std::string_view fileRoot, double mutRate,
                                   unsigned int samples) {
  if (!CSFSFile.empty() && fs::exists(CSFSFile)) {
    fmt::print("Precomputed CSFS will be loaded from file: {}\n", CSFSFile);
    auto csfs = CSFS::loadFromFile(CSFSFile);
    return calculateDecodingQuantities(csfs, demo, disc, fileRoot, freq, mutRate, samples);
  } else if (CSFSFile.empty()) {
    fmt::print("New CSFS will be calculated\n");
    return calculateCsfsAndDecodingQuantities(demo, disc, fileRoot, freq, mutRate, samples);
  } else {
    throw std::runtime_error(fmt::format("Specified CSFS file ({}) is invalid\n", CSFSFile));
  }
}

DecodingQuantities calculateCsfsAndDecodingQuantities(const Demography& demo, const Discretization& disc,
                                                      std::string_view fileRoot, const Frequencies& freq,
                                                      const double mutRate, const unsigned int samples) {

  // Get the array times and sizes, and remove the additional element added to the end of each array
  auto [arrayTime, arraySize] = getDemographicInfo(demo);
  arrayTime.pop_back();
  arraySize.pop_back();

  std::vector<double> arrayDisc = getDiscretizationInfo(disc, arrayTime, arraySize);

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

  return calculateDecodingQuantities(csfs, demo, disc, fileRoot, freq, mutRate, samples, std::move(arrayDisc));
}

std::tuple<std::vector<double>, std::vector<double>> getDemographicInfo(const Demography& demo) {
  std::vector<double> times;
  std::vector<double> sizes;

  if (demo.isFile()) {
    auto ts = readDemographic(demo.getDemography());
    times = ts.first;
    sizes = ts.second;
  } else {
    auto [t, s] = demo::getBuiltInDemography(demo.getDemography());
    times = std::move(t);
    sizes = std::move(s);
  }
  if (times.empty() || sizes.empty()) {
    throw std::runtime_error(fmt::format("Unknown error getting demography: {}", demo.getDemography()));
  }
  times.emplace_back(std::numeric_limits<double>::infinity());
  sizes.emplace_back(sizes.back()); // duplicate last element
  return std::make_tuple(times, sizes);
}

std::vector<double> getDiscretizationInfo(const Discretization& disc, const std::vector<double>& times,
                                          const std::vector<double>& sizes) {
  std::vector<double> discs;

  if (disc.isFile()) {
    fmt::print("Will read discretization intervals from {} ...\n", disc.getDiscretizationFile());
    discs = readDiscretization(disc.getDiscretizationFile());
  } else {
    // The pre-specified discretization quantiles
    discs = disc.getDiscretizationPoints();

    // If none pre-specified we can generate them all
    if (discs.empty()) {
      discs.emplace_back(0.0);
      fmt::print("Calculating {} discretization intervals from coalescent distribution.\n",
                 disc.getNumAdditionalPoints());
    } else {
      fmt::print("Using the following pre-specified discretization intervals: {}\n and calculating {} additional "
                 "intervals from coalescent distribution.\n",
                 discs, disc.getNumAdditionalPoints());
    }

    // We now shift by the final pre-specified quantile, and generate additional quantiles as required
    double lastPoint = discs.back();

    std::vector<double> newTimes;
    newTimes.reserve(times.size());
    std::vector<double> newSizes = sizes;

    for (double time : times) {
      newTimes.push_back(time - lastPoint);
    }

    if (newTimes.front() > 0.0 || newTimes.back() <= 0.0) {
      throw std::runtime_error("Something unexpected went wrong while calculating the discretization...\n");
    }

    auto i2 =
        std::adjacent_find(newTimes.begin(), newTimes.end(), [](double a, double b) { return a <= 0.0 && b > 0; });

    auto n = std::distance(newTimes.begin(), i2);

    newTimes.erase(newTimes.begin(), newTimes.begin() + n);
    newSizes.erase(newSizes.begin(), newSizes.begin() + n);
    newTimes.front() = 0.0;

    std::vector<double> newDiscs =
        Transition::getTimeExponentialQuantiles(1 + disc.getNumAdditionalPoints(), newTimes, newSizes);

    for (auto discNum = 1ul; discNum < newDiscs.size(); ++discNum) {
      discs.push_back(lastPoint + newDiscs.at(discNum));
    }
  }

  discs.emplace_back(std::numeric_limits<double>::infinity());

  return discs;
}

} // namespace asmc
