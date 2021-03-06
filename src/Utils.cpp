// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include "Utils.hpp"
#include "EigenTypes.hpp"
#include <array>
#include <cassert>
#include <cmath>
#include <string_view>
#include <sstream>
#include <fstream>

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

  char* successful_read = nullptr;
  do {
    successful_read = gzgets(gzFileHandle, buffer.data(), static_cast<int>(buffer.size()));

    if (successful_read != Z_NULL) {
      line += buffer.data();
    }

  } while (!line.empty() && line.back() != '\n' && !gzeof(gzFileHandle));

  if (!line.empty() && line.back() == '\n') {
    line.pop_back();
  }

  return line;
}

int writegz(gzFile& file, const std::string& s) {
  return gzwrite(file, s.c_str(), static_cast<unsigned int>(s.size()));
}

std::pair<std::vector<double>, std::vector<double>> readDemographic(std::string_view demographicFile) {
  std::vector<double> times, sizes;
  std::ifstream file(demographicFile.data());
  std::string line;
  double size = {}, t = {};
  while(std::getline(file, line)) {
    std::stringstream ss(line);
    ss >> t >> size;
    times.emplace_back(t);
    sizes.emplace_back(size);
  }
  return std::make_pair(times, sizes);
}

std::vector<double> readDiscretization(std::string_view discretizationFile) {
  std::vector<double> discs;
  std::ifstream file(discretizationFile.data());
  std::string line;
  while(std::getline(file, line)) discs.emplace_back(std::stod(line));
  return discs;
}

std::string matToString(const mat_dt& m, int precision) {
  std::stringstream ss;
  ss.precision(precision);
  if (!m.size()) return "";
  auto nrows = m.rows();
  auto ncols = m.cols();
  for(unsigned i = 0; i < nrows; i++) {
    for(unsigned j = 0; j < ncols - 1; j++) ss << m(i, j) << "\t";
    ss << m(i, ncols - 1) << "\n";
  }
  return ss.str();
}

std::string vecToString(const vec_dt& v, int precision) {
  std::stringstream ss;
  ss.precision(precision);
  if (!v.size()) return "";
  for(unsigned i = 0; i < v.size() - 1; i++) ss << v[i] << "\t";
  ss << v[v.size() - 1];
  return ss.str();
}

} // namespace asmc
