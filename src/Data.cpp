// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include "Data.hpp"

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <array>
#include <exception>
#include <filesystem>
#include <fstream>

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
  std::string line;
  std::stringstream tokens;

  // Fields in frequency file
  int chr, popSize;
  std::string SNP;
  char A1, A2;
  double freq;

  std::ifstream file;
  file.open(freqFile.data());
  if (file.is_open()) {
    getline(file, line);  // Ignore header
    while(getline(file, line)) {
      tokens = std::stringstream(line);
      tokens >> chr >> SNP >> A1 >> A2 >> freq >> popSize;
      if (popSize > mHaploidSampleSize)
        mHaploidSampleSize = popSize;
      mAllSNPsFreq << freq;
      mAllSNPsMinorAlleles << int(popSize * freq);
      mAllSNPsAlleleCounts << popSize;
    }
    file.close();
  } else {
    throw std::runtime_error(fmt::format("Could not read freq file: {}", freqFile));
  }
  fmt::print("{}", freqFile);
}

void Data::computeMinorAlleleFrequenciesFromHaps(std::string_view hapsFileRoot) {
  std::string hapsFileName = identifyAppropriateHapsFile(hapsFileRoot);
  std::stringstream tokens;
  std::string line, token;

  // TODO: Support .gz type
  int DAcount = 0;
  int samples = 0;
  double DAFreq;
  int pos = 0;
  int monomorphic = 0;

  std::ifstream hapsFile;
  hapsFile.open(hapsFileName.data());
  if(hapsFile.is_open()) {
    while (getline(hapsFile, line)) {
      tokens = std::stringstream(line);
      for(int i = 0; i < 5; i++)
        tokens >> token;  // ignore the first five tokens
      while (tokens >> token) {
        if (token == "1") DAcount++;
        samples++;
      }
      if (samples > mHaploidSampleSize)
        mHaploidSampleSize = samples;
      if (samples % 2 != 0)
        throw std::runtime_error("Haps file contains a line with odd haploid sample size.");
      if (DAcount > samples / 2)
        DAcount = samples - DAcount;
      DAFreq = DAcount / double(samples);
      mAllSNPsFreq << std::min(DAFreq, 1 - DAFreq);
      mAllSNPsMinorAlleles << DAcount;
      mAllSNPsAlleleCounts << samples;
      if (DAcount == 0)
        monomorphic++;
      pos++;
    }
    hapsFile.close();
    fmt::print("Computed frequencies for {} haploid samples and {} "
               "markers, {} of which are monomorphic.",
               mHaploidSampleSize, pos, monomorphic);

  } else {
    throw std::runtime_error(fmt::format("Could not open haps file: {}", hapsFileName));
  }
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
