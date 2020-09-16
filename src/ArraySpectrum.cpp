// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include <vector>
#include <numeric>
#include <iostream>
#include "Data.hpp"
#include "Utils.hpp"
#include "ArraySpectrum.hpp"

namespace asmc {

ArraySpectrum::ArraySpectrum(Data data, unsigned samples) {
  std::unordered_map<double, std::tuple<int, int, int>> distributions;
  std::unordered_map<double, int> distCounts;
  std::vector<double> spectrum;

  // get hypergeometric for each allele
  int monoMorphic = 0;
  auto alleleCounts = data.getAllSNPsAlleleCounts();
  auto minorAlleles = data.getAllSNPsMinorAlleles();
  auto freqs = data.getAllSNPsFreq();
  for (int i = 0; i < alleleCounts.size(); i++) {
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
  for (std::pair<double, int> distCount: distCounts) {
    for (unsigned i = 0; i <= samples; i++)
      spectrum.push_back(
          std::apply(hypergeometricPmf, std::tuple_cat(distributions[distCount.first], std::tie(i)))
          * distCount.second);
  }
  // the term in 0 will contain monomorphic samples that are either present in the data, or due to subsampling
  spectrum[0] += monoMorphic;

  normalize(spectrum);

  // store monomorphic probability separately, then renormalize excluding monomorphic at 0 and samples
  // add alleles for which all samples are carriers to monomorphic probability
  mMonomorphicProbability = spectrum[0] + spectrum[samples];
  std::cout << "Probability of a site being monomorphic due to subsampling: "
    << std::round(mMonomorphicProbability * 1000) / 1000.0 << std::endl;
  // renormalize without monomorphic
  spectrum[0] = 0.;
  spectrum[samples] = 0.;
  normalize(spectrum);

  // fold to minor allele
  unsigned halfTotal = samples / 2;
  for (unsigned i = 0; i < halfTotal; i++) {
      spectrum[i] = spectrum[i] + spectrum[samples - i];
  }
  mSpectrum = std::vector<double>(spectrum.begin(), spectrum.begin() + halfTotal + 1);
}

void normalize(std::vector<double>& spectrum) {
  double tot = std::reduce(spectrum.begin(), spectrum.end());
  std::for_each(spectrum.begin(), spectrum.end(), [tot](double& s) { s /= tot;});
}

} // namespace asmc

