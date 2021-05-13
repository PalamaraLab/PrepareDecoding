// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include "DecodingQuantities.hpp"
#include "PrepareDecoding.hpp"
#include "ThinParameterTypes.hpp"
#include "Utils.hpp"
#include <cstdlib>
#include <cxxopts.hpp>
#include <filesystem>
#include <fmt/core.h>

namespace fs = std::filesystem;
using namespace asmc;

/* DecodingQuantities prepareDecoding(std::string_view demographicFile, std::string_view discretizationFile, */
/*                                           int coalescentQuantiles, int mutationAgeIntervals, std::string_view
 * fileRoot, */
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

  cxxopts::Options options("ASMCprepareDecoding", "Precompute decoding quantities for ASMC");

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
       cxxopts::value<std::string>()->default_value(""))
      ("q,coalescentQuantiles", "Desired number of discretization intervals "
       " (quantiles from the pairwise coalescent distribution)",
       cxxopts::value<int>()->default_value("0"))
      ("o,outputFileRoot", "Output file root.",
       cxxopts::value<std::string>()->default_value(""))
      ("f,fileRoot", "Root for name of hap/samples files from which allele "
                     " frequencies for array data should be read.",
       cxxopts::value<std::string>()->default_value(""))
      ("F,freqFile", "Plink .frq file containing minor allele frequencies for array data.",
       cxxopts::value<std::string>()->default_value(""))
  ;
  // clang-format on

  auto result = options.parse(argc, argv);

  const std::string demographicFile = result["demography"].as<std::string>();
  const std::string discretizationFile = result["discretization"].as<std::string>();
  const std::string CSFSFile = result["CSFS"].as<std::string>();
  const std::string freqFile = result["freqFile"].as<std::string>();
  const std::string fileRoot = result["fileRoot"].as<std::string>();
  const std::string outputFileRoot = result["outputFileRoot"].as<std::string>();

  const int coalescentQuantiles = result["coalescentQuantiles"].as<int>();
  const int mutationAgeIntervals = result["mutationQuantiles"].as<int>();

  const unsigned int samples = result["samples"].as<unsigned int>();

  const double mutRate = result["mut"].as<double>();

  DecodingQuantities dq =
      CSFSFile.empty()
          ? calculateCsfsAndPrepareDecoding(Demography(demographicFile), discretizationFile, coalescentQuantiles,
                                            mutationAgeIntervals, fileRoot, freqFile, mutRate, samples)
          : prepareDecodingPrecalculatedCsfs(CSFSFile, Demography(demographicFile), discretizationFile,
                                             coalescentQuantiles, mutationAgeIntervals, fileRoot, freqFile, mutRate,
                                             samples);

  dq.saveDecodingQuantities(outputFileRoot);
  dq.saveIntervals(outputFileRoot);
  dq.saveCsfs(outputFileRoot);
  fmt::print("Done\n");
}


