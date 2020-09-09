// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include <cxxopts.hpp>
#include <fmt/core.h>

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
       cxxopts::value<int>()->default_value("0"))
      ("m,mut", "Mutation rate assumed by demographic model",
       cxxopts::value<double>()->default_value("1.65e-8"))
      ("d,discretization", "File with vector of time discretization intervals",
       cxxopts::value<std::string>())
      ("q,coalescentQuantiles", "Desired number of discretization intervals "
       " (quantiles from the pairwise coalescent distribution)",
       cxxopts::value<std::string>())
      ("o,OutputFileRoot", "Output file root.",
       cxxopts::value<std::string>())
      ("f,fileRoot", "Root for name of hap/samples files from which allele "
                     " frequencies for array data should be read.",
       cxxopts::value<std::string>())
      ("F,freqFile", "Plink .frq file containing minor allele frequencies for array data.",
       cxxopts::value<std::string>())
      ;
  options.parse(argc, argv);
}
