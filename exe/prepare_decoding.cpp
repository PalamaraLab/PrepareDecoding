// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include <cxxopts.hpp>
#include <fmt/core.h>
#include <filesystem>
#include <cstdlib>
#include "Transition.hpp"
#include "Utils.hpp"
#include "Data.hpp"
#include "Csfs.hpp"
#include "DecodingQuantities.hpp"

namespace fs = std::filesystem;
using namespace asmc;

int main(int argc, char* argv[]) {

  cxxopts::Options options("ASMCprepareDecoding",
          "Precompute decoding quantities for ASMC");

  // clang-format off
  options.add_options()
      ("D,demography",
       "File with demographic model. If not specified, the default EU model will be used.",
       cxxopts::value<std::string>()->default_value(""))
      ("C,CSFS",
       "File with precomputed CSFS. If not specified this program will try to compute it internally.",
       cxxopts::value<std::string>()->default_value(""))
      ("n,samples",
       "Number of total samples to be used for emission probability. If not specified "
       "will use min(300, sample size)",
       cxxopts::value<unsigned int>()->default_value("300"))
      ("m,mut", "Mutation rate assumed by demographic model",
       cxxopts::value<double>()->default_value("1.65e-8"))
      ("mutationQuantiles", "Infer desired number of discretization intervals from the distribution of mutation age.",
       cxxopts::value<int>()->default_value("0"))
      ("d,discretization", "File with vector of time discretization intervals",
       cxxopts::value<std::string>())
      ("q,coalescentQuantiles", "Desired number of discretization intervals "
       " (quantiles from the pairwise coalescent distribution)",
       cxxopts::value<int>())
      ("o,OutputFileRoot", "Output file root.",
       cxxopts::value<std::string>())
      ("f,fileRoot", "Root for name of hap/samples files from which allele "
                     " frequencies for array data should be read.",
       cxxopts::value<std::string>())
      ("F,freqFile", "Plink .frq file containing minor allele frequencies for array data.",
       cxxopts::value<std::string>())
      ;
  auto result = options.parse(argc, argv);

  std::string demographicFile;
  std::vector<double> times, sizes;

  // Read emographic model
  if(result.count("demography")) demographicFile = result["demography"].as<std::string>();
  if(!demographicFile.empty() & fs::exists(demographicFile)) {
    fmt::print("Will read demographic model from {} ...\n", demographicFile);
    auto ts = readDemographic(demographicFile);
    times = ts.first;
    sizes = ts.second;
  } else {
    times = Transition::EUtime;
    sizes = Transition::EUsize;
    fmt::print("Did not input a demographic model, using default EU model.\n");
  }
  times.emplace_back(std::numeric_limits<double>::infinity());
  sizes.emplace_back(sizes.back()); // duplicate last element

  std::vector<double> discs;
  std::string discretizationFile;
  int coalescentQuantiles = 0, mutationAgeIntervals = 0;
  if(result.count("discretization")) discretizationFile = result["discretization"].as<std::string>();
  if(result.count("coalescentQuantiles")) coalescentQuantiles = result["coalescentQuantiles"].as<int>();
  if(result.count("mutationQuantiles")) mutationAgeIntervals = result["mutationQuantiles"].as<int>();
  
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

  std::string fileRoot, freqFile;
  Data data;
  if(result.count("fileRoot")) fileRoot = result["fileRoot"].as<std::string>();
  if(result.count("freqFile")) freqFile = result["freqFile"].as<std::string>();
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
  auto mutRate = result["mut"].as<double>();
  fmt::print("Will use mutation rate mu = {}.", mutRate);
  auto samples = result["samples"].as<unsigned int>();
  samples = std::min(samples, data.getHaploidSampleSize());
  fmt::print("Number of samples in CSFS calculations: {}.\n", samples);

  // Transition
  Transition transition(times, sizes, discs, CSC);

  // Parse CSFS
  std::string CSFSFile;
  if(result.count("CSFS")) CSFSFile = result["CSFS"].as<std::string>();
  CSFS csfs;
  if(!CSFSFile.empty() & fs::exists(CSFSFile)) {
    fmt::print("Will load precomputed CSFS from {} ...\n", CSFSFile);
    csfs = CSFS::loadFromFile(CSFSFile);
    fmt::print("Verifying CSFS loaded from {} ...", CSFSFile);
    if(csfs.verify(times, sizes, mutRate, samples, discs)) {
      fmt::print("Verified " + std::to_string(csfs.getCSFS().size()) + " CSFS entries.\n");
    } else {
      fmt::print(
          "CSFS could not be verified, aborting.\n"
          "You can use the Python version of PrepareDecoding to fetch CSFS\n"
          "and prepare decoding quantities.\n");
      return EXIT_FAILURE;
    }
  } else {
    throw std::runtime_error("Valid CSFS file needs to be specified\n");
  }
  csfs.fixAscertainment(data, samples, transition);

  if(!result.count("outputFileRoot")) throw std::runtime_error("Required option --outputFileRoot missing.\n");
  auto outputFileRoot = result["outputFileRoot"].as<std::string>();

  // Build decoding quantities
  fmt::print("\nBuilding decoding quantities...\n");
  DecodingQuantities decodingQuantities(csfs, transition, mutRate);
  decodingQuantities.saveDecodingQuantities(outputFileRoot);
  decodingQuantities.saveIntervals(outputFileRoot);
  fmt::print("Done\n");

}


