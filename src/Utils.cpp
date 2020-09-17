// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include "Utils.hpp"

#include <array>
#include <cassert>
#include <cmath>

namespace asmc {

double hypergeometricPmf(const int populationSize, const int numberOfSuccesses, const int sampleSize,
                         const int observedSuccesses) {

  assert(populationSize >= 0);
  assert(numberOfSuccesses >= 0);
  assert(sampleSize >= 0);
  assert(observedSuccesses >= 0);

  assert(populationSize >= numberOfSuccesses);

  assert(observedSuccesses <= sampleSize);
  assert(observedSuccesses <= numberOfSuccesses);

  const long double num1 = std::lgamma(1.L + numberOfSuccesses);
  const long double num2 = std::lgamma(1.L + sampleSize);
  const long double num3 = std::lgamma(1.L + populationSize - numberOfSuccesses);
  const long double num4 = std::lgamma(1.L + populationSize - sampleSize);

  const long double den1 = std::lgamma(1.L + observedSuccesses);
  const long double den2 = std::lgamma(1.L + populationSize);
  const long double den3 = std::lgamma(1.L + numberOfSuccesses - observedSuccesses);
  const long double den4 = std::lgamma(1.L + sampleSize - observedSuccesses);
  const long double den5 = std::lgamma(1.L + populationSize + observedSuccesses - numberOfSuccesses - sampleSize);

  return static_cast<double>(std::exp(num1 + num2 + num3 + num4 - den1 - den2 - den3 - den4 - den5));
}

std::string readNextLineFromGzip(gzFile& gzFileHandle) {

  std::array<char, 512> buffer = {};
  std::string line;

  {
    char* successful_read = nullptr;
    do {
      successful_read = gzgets(gzFileHandle, buffer.data(), buffer.size());

      if (successful_read != Z_NULL) {
        line += buffer.data();
      }

    } while (line.back() != '\n' && !gzeof(gzFileHandle));
  }

  if(line.back() == '\n') {
    line.pop_back();
  }

  return line;
}

} // namespace asmc
