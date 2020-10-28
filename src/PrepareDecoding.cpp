// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

// Main PrepareDecoding functions

#include <fmt/core.h>
#include <filesystem>
#include <cstdlib>
#include "Transition.hpp"
#include "Utils.hpp"
#include "Data.hpp"
#include "Csfs.hpp"
#include "DecodingQuantities.hpp"
#include "PrepareDecoding.hpp"

namespace fs = std::filesystem;

namespace asmc {

DecodingQuantities prepareDecoding(std::string_view CSFSFile, std::string_view demographicFile,
                                   std::string_view discretizationFile, int coalescentQuantiles,
                                   int mutationAgeIntervals, std::string_view fileRoot, std::string_view freqFile,
                                   double mutRate, unsigned int samples) {

  std::vector<double> times, sizes, discs;
  Data data;
  // Read emographic model
  if(!demographicFile.empty() & fs::exists(demographicFile)) {
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

  // Read discretization
  if(!discretizationFile.empty() & fs::exists(discretizationFile)) {
    fmt::print("Will read discretization intervals from {} ...\n", discretizationFile);
    discs = readDiscretization(discretizationFile);
  } else if (coalescentQuantiles > 0) {
    discs = Transition::getTimeExponentialQuantiles(coalescentQuantiles, times, sizes);
    fmt::print("Using {} discretization intervals from coalescent distribution.\n", coalescentQuantiles);
  } else if (mutationAgeIntervals > 0) {
    discs = Transition::getTimeErlangQuantiles(mutationAgeIntervals, times, sizes);
    fmt::print("Using {} discretization intervals from coalescent distribution.\n", mutationAgeIntervals);
  }
  if(discs.empty())
    throw std::runtime_error("Specify a valid option from --discretization, --coalescentQuantiles, --mutationQuantiles\n");
  discs.emplace_back(std::numeric_limits<double>::infinity());

  if(!fileRoot.empty()) fmt::print("Files will be read from: {}*\n", fileRoot);
  if(!freqFile.empty()) {
    if(!fs::exists(freqFile))
      throw std::runtime_error(fmt::format("Could not open {}\n", freqFile));
    fmt::print("Will load minor allele frequencies from {} ...\n", freqFile);
    data.addFreq(freqFile);
  } else {
    if(fileRoot.empty()) throw std::runtime_error("Either one of --freqFile or --fileRoot has to be specified\n");
    fmt::print("Will compute allele frequencies in haps files with root {}.", fileRoot);
    data = Data(fileRoot);
  }
  fmt::print("Will use mutation rate mu = {}.\n", mutRate);
  samples = std::min(samples, data.getHaploidSampleSize());
  fmt::print("Number of samples in CSFS calculations: {}.\n", samples);

  // Transition
  Transition transition(times, sizes, discs, CSC);

  // Parse CSFS
  CSFS csfs;
  if(!CSFSFile.empty() & fs::exists(CSFSFile)) {
    fmt::print("Will load precomputed CSFS from {} ...\n", CSFSFile);
    csfs = CSFS::loadFromFile(CSFSFile);
    fmt::print("Verifying CSFS loaded from {} ...\n", CSFSFile);
    if(csfs.verify(times, sizes, mutRate, samples, discs)) {
      fmt::print("Verified " + std::to_string(csfs.getCSFS().size()) + " CSFS entries.\n");
    } else {
      throw std::runtime_error("CSFS could not be verified. The Python version can be used"
          "to fetch CSFS and prepare decoding quantities.");
    }
  } else {
    throw std::runtime_error("Valid CSFS file needs to be specified\n");
  }
  csfs.fixAscertainment(data, samples, transition);

  // Build decoding quantities
  fmt::print("\nBuilding decoding quantities...\n");
  return DecodingQuantities(csfs, transition, mutRate);
}

} // namespace asmc
