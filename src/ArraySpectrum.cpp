// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include <vector>
#include <iostream>
#include <unordered_map>
#include "Data.hpp"
#include "Utils.hpp"
#include "ArraySpectrum.hpp"
#include "DefaultArraySpectra.hpp"

#include <fmt/format.h>
#include <fmt/ostream.h>

namespace asmc {

ArraySpectrum::ArraySpectrum(Data data, unsigned samples) {

  // If we're using a built-in, simply generate from the existing values
  if(const Frequencies& freq = data.getFreq(); freq.isBuiltIn()) {
    assert(freq.getNumSamples() == samples);

    std::vector<double> spectrumCopy;
    if(freq.getFreqIdentifier() == "UKBB") {
      switch (freq.getNumSamples()) {
      case 50:
        spectrumCopy = std::vector<double>(spectra::spectrumUkbb50.begin(), spectra::spectrumUkbb50.end());
        mMonomorphicProbability = spectra::pssUkbb50;
        break;
      case 100:
        spectrumCopy = std::vector<double>(spectra::spectrumUkbb100.begin(), spectra::spectrumUkbb100.end());
        mMonomorphicProbability = spectra::pssUkbb100;
        break;
      case 200:
        spectrumCopy = std::vector<double>(spectra::spectrumUkbb200.begin(), spectra::spectrumUkbb200.end());
        mMonomorphicProbability = spectra::pssUkbb200;
        break;
      case 300:
        spectrumCopy = std::vector<double>(spectra::spectrumUkbb300.begin(), spectra::spectrumUkbb300.end());
        mMonomorphicProbability = spectra::pssUkbb300;
        break;
      }
    }

    mSpectrum.resize(static_cast<Eigen::Index>(spectrumCopy.size()));
    for (auto i = 0ul; i < spectrumCopy.size(); ++i) {
      mSpectrum[static_cast<Eigen::Index>(i)] = spectrumCopy.at(i);
    }

    return;
  }

  // If not using a built-in, we have to calculate them
  std::unordered_map<double, std::tuple<int, int, int>> distributions;
  std::unordered_map<double, int> distCounts;
  array_dt spectrum(samples + 1);
  spectrum.setZero();

  // get hypergeometric for each allele
  int monoMorphic = 0;
  auto alleleCounts = data.getAllSNPsAlleleCounts();
  auto minorAlleles = data.getAllSNPsMinorAlleles();
  auto freqs = data.getAllSNPsFreq();
  for (unsigned i = 0; i < alleleCounts.size(); i++) {
      double freq = freqs[i];
      if (minorAlleles[i] == 0) {
          monoMorphic++;
          continue;
      }
      auto found = distributions.find(freq);
      if (found == distributions.end()) {
        distributions.emplace(freq,
            std::make_tuple(alleleCounts[i], minorAlleles[i], samples));
        distCounts.emplace(freq, 1);
      } else {
        distCounts[freq]++;
      }
  }

  // sum all hypergeometrics to get spectrum in subsample (exclude [0] and [samples], which are monomorphic.
  for (auto [_freq, _count] : distCounts) {
    for (unsigned i = 0; i <= samples; i++) {
      auto [p0, p1, p2] = distributions[_freq];
      spectrum[i] += asmc::hypergeometricPmf(p0, p1, p2, static_cast<int>(i)) * _count;
    }
  }
  // the term in 0 will contain monomorphic samples that are either present in the data, or due to subsampling
  spectrum[0] += monoMorphic;
  spectrum /= spectrum.sum();

  // store monomorphic probability separately, then renormalize excluding monomorphic at 0 and samples
  // add alleles for which all samples are carriers to monomorphic probability
  mMonomorphicProbability = spectrum[0] + spectrum[samples];
  fmt::print("Probability of a site being monomorphic due to subsampling: {:.15f}\n", mMonomorphicProbability);

  // renormalize without monomorphic
  spectrum[0] = 0.;
  spectrum[samples] = 0.;
  spectrum /= spectrum.sum();

  // fold to minor allele
  unsigned halfTotal = samples / 2;
  for (unsigned i = 0; i < halfTotal; i++) spectrum[i] += spectrum[samples - i];
  mSpectrum = spectrum.head(halfTotal + 1);
}

} // namespace asmc

