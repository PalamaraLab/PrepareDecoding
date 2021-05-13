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
#include <cstdlib>
#include <filesystem>
#include <fmt/core.h>
#include <utility>

#include "smcpp.hpp"

namespace fs = std::filesystem;

namespace asmc {

DecodingQuantities prepareDecoding(CSFS& csfs, const Demography& demo, std::string_view discretizationFile,
                                   int coalescentQuantiles, int mutationAgeIntervals, std::string_view fileRoot,
                                   std::string_view freqFile, double mutRate, unsigned int samples) {

  auto [times, sizes] = getDemographicInfo(demo);
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

DecodingQuantities calculateCsfsAndPrepareDecoding(const Demography& demo, std::string_view discretizationFile,
                                                   const int coalescentQuantiles, const int mutationAgeIntervals,
                                                   std::string_view fileRoot, std::string_view freqFile,
                                                   const double mutRate, const unsigned int samples) {

  // Get the array times and sizes, and remove the additional element added to the end of each array
  auto [arrayTime, arraySize] = getDemographicInfo(demo);
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

  return prepareDecoding(csfs, demo, discretizationFile, coalescentQuantiles, mutationAgeIntervals, fileRoot,
                         freqFile, mutRate, samples);
}

DecodingQuantities prepareDecodingPrecalculatedCsfs(std::string_view CSFSFile, const Demography& demo,
                                                    std::string_view discretizationFile, int coalescentQuantiles,
                                                    int mutationAgeIntervals, std::string_view fileRoot,
                                                    std::string_view freqFile, double mutRate, unsigned int samples) {
  if (!CSFSFile.empty() & fs::exists(CSFSFile)) {
    fmt::print("Will load precomputed CSFS from {} ...\n", CSFSFile);
    auto csfs = CSFS::loadFromFile(CSFSFile);
    return prepareDecoding(csfs, demo, discretizationFile, coalescentQuantiles, mutationAgeIntervals,
                           fileRoot, freqFile, mutRate, samples);
  } else {
    throw std::runtime_error("Valid CSFS file needs to be specified\n");
  }
}

std::tuple<std::vector<double>, std::vector<double>> getDemographicInfo(const Demography& demo) {
  std::vector<double> times;
  std::vector<double> sizes;

  if (demo.isFile()) {
    auto ts = readDemographic(demo.getDemography());
    times = ts.first;
    sizes = ts.second;
  } else {
    if (demo.getDemography() == "ACB") {
      times = std::vector<double>(demo::timesACB.begin(), demo::timesACB.end());
      sizes = std::vector<double>(demo::sizesACB.begin(), demo::sizesACB.end());
    } else if (demo.getDemography() == "ASW") {
      times = std::vector<double>(demo::timesASW.begin(), demo::timesASW.end());
      sizes = std::vector<double>(demo::sizesASW.begin(), demo::sizesASW.end());
    } else if (demo.getDemography() == "BEB") {
      times = std::vector<double>(demo::timesBEB.begin(), demo::timesBEB.end());
      sizes = std::vector<double>(demo::sizesBEB.begin(), demo::sizesBEB.end());
    } else if (demo.getDemography() == "CDX") {
      times = std::vector<double>(demo::timesCDX.begin(), demo::timesCDX.end());
      sizes = std::vector<double>(demo::sizesCDX.begin(), demo::sizesCDX.end());
    } else if (demo.getDemography() == "CEU") {
      times = std::vector<double>(demo::timesCEU.begin(), demo::timesCEU.end());
      sizes = std::vector<double>(demo::sizesCEU.begin(), demo::sizesCEU.end());
    } else if (demo.getDemography() == "CHB") {
      times = std::vector<double>(demo::timesCHB.begin(), demo::timesCHB.end());
      sizes = std::vector<double>(demo::sizesCHB.begin(), demo::sizesCHB.end());
    } else if (demo.getDemography() == "CHS") {
      times = std::vector<double>(demo::timesCHS.begin(), demo::timesCHS.end());
      sizes = std::vector<double>(demo::sizesCHS.begin(), demo::sizesCHS.end());
    } else if (demo.getDemography() == "CLM") {
      times = std::vector<double>(demo::timesCLM.begin(), demo::timesCLM.end());
      sizes = std::vector<double>(demo::sizesCLM.begin(), demo::sizesCLM.end());
    } else if (demo.getDemography() == "ESN") {
      times = std::vector<double>(demo::timesESN.begin(), demo::timesESN.end());
      sizes = std::vector<double>(demo::sizesESN.begin(), demo::sizesESN.end());
    } else if (demo.getDemography() == "FIN") {
      times = std::vector<double>(demo::timesFIN.begin(), demo::timesFIN.end());
      sizes = std::vector<double>(demo::sizesFIN.begin(), demo::sizesFIN.end());
    } else if (demo.getDemography() == "GBR") {
      times = std::vector<double>(demo::timesGBR.begin(), demo::timesGBR.end());
      sizes = std::vector<double>(demo::sizesGBR.begin(), demo::sizesGBR.end());
    } else if (demo.getDemography() == "GIH") {
      times = std::vector<double>(demo::timesGIH.begin(), demo::timesGIH.end());
      sizes = std::vector<double>(demo::sizesGIH.begin(), demo::sizesGIH.end());
    } else if (demo.getDemography() == "GWD") {
      times = std::vector<double>(demo::timesGWD.begin(), demo::timesGWD.end());
      sizes = std::vector<double>(demo::sizesGWD.begin(), demo::sizesGWD.end());
    } else if (demo.getDemography() == "IBS") {
      times = std::vector<double>(demo::timesIBS.begin(), demo::timesIBS.end());
      sizes = std::vector<double>(demo::sizesIBS.begin(), demo::sizesIBS.end());
    } else if (demo.getDemography() == "ITU") {
      times = std::vector<double>(demo::timesITU.begin(), demo::timesITU.end());
      sizes = std::vector<double>(demo::sizesITU.begin(), demo::sizesITU.end());
    } else if (demo.getDemography() == "JPT") {
      times = std::vector<double>(demo::timesJPT.begin(), demo::timesJPT.end());
      sizes = std::vector<double>(demo::sizesJPT.begin(), demo::sizesJPT.end());
    } else if (demo.getDemography() == "KHV") {
      times = std::vector<double>(demo::timesKHV.begin(), demo::timesKHV.end());
      sizes = std::vector<double>(demo::sizesKHV.begin(), demo::sizesKHV.end());
    } else if (demo.getDemography() == "LWK") {
      times = std::vector<double>(demo::timesLWK.begin(), demo::timesLWK.end());
      sizes = std::vector<double>(demo::sizesLWK.begin(), demo::sizesLWK.end());
    } else if (demo.getDemography() == "MSL") {
      times = std::vector<double>(demo::timesMSL.begin(), demo::timesMSL.end());
      sizes = std::vector<double>(demo::sizesMSL.begin(), demo::sizesMSL.end());
    } else if (demo.getDemography() == "MXL") {
      times = std::vector<double>(demo::timesMXL.begin(), demo::timesMXL.end());
      sizes = std::vector<double>(demo::sizesMXL.begin(), demo::sizesMXL.end());
    } else if (demo.getDemography() == "PEL") {
      times = std::vector<double>(demo::timesPEL.begin(), demo::timesPEL.end());
      sizes = std::vector<double>(demo::sizesPEL.begin(), demo::sizesPEL.end());
    } else if (demo.getDemography() == "PJL") {
      times = std::vector<double>(demo::timesPJL.begin(), demo::timesPJL.end());
      sizes = std::vector<double>(demo::sizesPJL.begin(), demo::sizesPJL.end());
    } else if (demo.getDemography() == "PUR") {
      times = std::vector<double>(demo::timesPUR.begin(), demo::timesPUR.end());
      sizes = std::vector<double>(demo::sizesPUR.begin(), demo::sizesPUR.end());
    } else if (demo.getDemography() == "STU") {
      times = std::vector<double>(demo::timesSTU.begin(), demo::timesSTU.end());
      sizes = std::vector<double>(demo::sizesSTU.begin(), demo::sizesSTU.end());
    } else if (demo.getDemography() == "TSI") {
      times = std::vector<double>(demo::timesTSI.begin(), demo::timesTSI.end());
      sizes = std::vector<double>(demo::sizesTSI.begin(), demo::sizesTSI.end());
    } else if (demo.getDemography() == "YRI") {
      times = std::vector<double>(demo::timesYRI.begin(), demo::timesYRI.end());
      sizes = std::vector<double>(demo::sizesYRI.begin(), demo::sizesYRI.end());
    }
  }
  if (times.empty() || sizes.empty()) {
    throw std::runtime_error(fmt::format("Unknown error getting demography: {}", demo.getDemography()));
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
