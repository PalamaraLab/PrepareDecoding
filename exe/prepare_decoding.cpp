// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include <cxxopts.hpp>
#include <fmt/core.h>
#include <filesystem>
#include <cstdlib>
#include "Utils.hpp"
#include "DecodingQuantities.hpp"
#include "PrepareDecoding.hpp"

namespace fs = std::filesystem;
using namespace asmc;

/* DecodingQuantities prepareDecoding(std::string_view demographicFile, std::string_view discretizationFile, */
/*                                           int coalescentQuantiles, int mutationAgeIntervals, std::string_view fileRoot, */
/*                                           std::string_view freqFile, double mutRate, unsigned int samples, */
/*                                           std::string_view CSFSFile) { */
int main(int argc, char* argv[]) {

  const std::string VERSION("1.0");
  const std::string VERSION_DATE("July 1, 2018");
  const std::string YEAR("2018");
  const std::string LICENSE("GPL v3");
  const std::string WEBSITE("www.palamaralab.org/software/ASMC");
  const std::string PROGRAM("ASMC");

  fmt::print(R"header(
 █████╗   ███████╗  ███╗   ███╗   ██████╗
██╔══██╗  ██╔════╝  ████╗ ████║  ██╔════╝
███████║  ███████╗  ██╔████╔██║  ██║
██╔══██║  ╚════██║  ██║╚██╔╝██║  ██║
██║  ██║  ███████║  ██║ ╚═╝ ██║  ╚██████╗
╚═╝  ╚═╝  ╚══════╝  ╚═╝     ╚═╝   ╚═════╝
)header");

  fmt::print("\nAscertained Sequentially Markovian Coalescent (ASMC) - Precompute decoding quantities v.{}\n", VERSION);
  fmt::print("GNU GPL v3, Copyright (C) {} Pier Palamara\n", YEAR);
  fmt::print("Manual: {}\n\n", WEBSITE);

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
      ("o,outputFileRoot", "Output file root.",
       cxxopts::value<std::string>())
      ("f,fileRoot", "Root for name of hap/samples files from which allele "
                     " frequencies for array data should be read.",
       cxxopts::value<std::string>())
      ("F,freqFile", "Plink .frq file containing minor allele frequencies for array data.",
       cxxopts::value<std::string>())
      ;
  auto result = options.parse(argc, argv);

  std::string demographicFile, discretizationFile, CSFSFile, fileRoot, freqFile;
  int coalescentQuantiles = 0;
  // Read emographic model
  if(result.count("demography")) demographicFile = result["demography"].as<std::string>();
  if(result.count("discretization")) discretizationFile = result["discretization"].as<std::string>();
  if(result.count("coalescentQuantiles")) coalescentQuantiles = result["coalescentQuantiles"].as<int>();
  auto mutationAgeIntervals = result["mutationQuantiles"].as<int>();
  if(result.count("fileRoot")) fileRoot = result["fileRoot"].as<std::string>();
  if(result.count("freqFile")) freqFile = result["freqFile"].as<std::string>();
  auto mutRate = result["mut"].as<double>();
  auto samples = result["samples"].as<unsigned int>();
  auto outputFileRoot = result["outputFileRoot"].as<std::string>();
  if(result.count("CSFS")) CSFSFile = result["CSFS"].as<std::string>();

  auto dq = prepareDecodingCSFSFile(CSFSFile, demographicFile, discretizationFile,
                            coalescentQuantiles, mutationAgeIntervals, fileRoot,
                            freqFile, mutRate, samples);
  dq.saveDecodingQuantities(outputFileRoot);
  dq.saveIntervals(outputFileRoot);
  fmt::print("Done\n");

}


