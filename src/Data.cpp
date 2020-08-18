// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include "Data.hpp"

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <array>
#include <exception>
#include <filesystem>

namespace asmc {

namespace fs = std::filesystem;

Data::Data(std::string_view hapsFileRoot) {

  if (auto frq_gz = fmt::format("{}.frq.gz", hapsFileRoot); fs::exists(frq_gz)) {
    readMinorAlleleFrequencies(frq_gz);
  } else if (auto frq = fmt::format("{}.frq", hapsFileRoot); fs::exists(frq)) {
    readMinorAlleleFrequencies(frq);
  } else {
    computeMinorAlleleFrequenciesFromHaps(hapsFileRoot);
  }
}

void Data::addFreq(std::string_view freqFile) {
  readMinorAlleleFrequencies(freqFile);
}

void Data::readMinorAlleleFrequencies(std::string_view freqFile) {
  fmt::print("{}", freqFile);
}

void Data::computeMinorAlleleFrequenciesFromHaps(std::string_view hapsFileRoot) {
  std::string hapsFile = identifyAppropriateHapsFile(hapsFileRoot);
}

std::string Data::identifyAppropriateHapsFile(std::string_view hapsFileRoot) {
  const std::array<std::string, 4> permittedHapsExt = {".hap.gz", ".hap", ".haps.gz", ".haps"};

  auto use_ext = std::find_if(permittedHapsExt.begin(), permittedHapsExt.end(), [hapsFileRoot](std::string_view ext) {
    return fs::exists(fmt::format("{}{}", hapsFileRoot, ext));
  });

  if (use_ext == permittedHapsExt.end()) {
    throw std::runtime_error(fmt::format("No haps file found at {} with any of the following extensions: {}",
                                         hapsFileRoot, permittedHapsExt));
  }

  return fmt::format("{}{}", hapsFileRoot, *use_ext);
}

} // namespace asmc
