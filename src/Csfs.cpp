// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include "Csfs.hpp"

#include "EigenTypes.hpp"

#include <iostream>
#include <fstream>
#include <fmt/core.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>

#include <stdexcept>

namespace asmc {


CSFS::CSFS(std::map<double, CSFSEntry> CSFS_) : mCSFS(std::move(CSFS_)),
  mAscertainedCSFS({}), mFoldedAscertainedCSFS({}) {
  if (!mCSFS.empty()) mFoldedCSFS = foldCSFS(mCSFS);
  stateMap.emplace(std::make_pair("Size:", CSFSParserState::Size));
  stateMap.emplace(std::make_pair("Time:", CSFSParserState::Time));
  stateMap.emplace(std::make_pair("Mu:", CSFSParserState::Mu));
  stateMap.emplace(std::make_pair("Samples:", CSFSParserState::Samples));
  stateMap.emplace(std::make_pair("Interval:", CSFSParserState::Interval));
}

CSFSParserState CSFS::currentState(const std::string& line) {
  std::stringstream lineStream(line);
  std::string token;
  lineStream >> token;
  auto search = stateMap.find(token);
  return ((search != stateMap.end()) ? search->second : CSFSParserState::CSFS);
}

std::pair<CSFSParserState, int> CSFS::nextState(CSFSParserState state, int line = 0) {
  switch (state) {
    case CSFSParserState::Null:
      return std::make_pair(CSFSParserState::Time, 0);

    case CSFSParserState::Time:
      return std::make_pair(CSFSParserState::Size, 0);

    case CSFSParserState::Size:
      return std::make_pair(CSFSParserState::Mu, 0);

    case CSFSParserState::Mu:
      return std::make_pair(CSFSParserState::Samples, 0);

    case CSFSParserState::Samples:
      return std::make_pair(CSFSParserState::Interval, 0);

    case CSFSParserState::Interval:
      return std::make_pair(CSFSParserState::CSFS, 0);

    case CSFSParserState::CSFS:
      if (line == 3) return std::make_pair(CSFSParserState::Time, 0);
      return std::make_pair(CSFSParserState::CSFS, line + 1);
  }
  return std::make_pair(CSFSParserState::Null, 0);
}

CSFS CSFS::loadFromFile(std::string_view filename) {
  std::ifstream file;
  file.open(filename.data());
  std::string line;
  std::map<double, CSFSEntry> parsed;
  std::vector<double> timeVector = {};
  std::vector<double> sizeVector = {};
  mat_dt csfs;
  std::string token;
  double t_dt = {}, from = {}, to = {}, mu = {};
  int samples = {};
  auto [state, subline] = nextState(CSFSParserState::Null);
  while(std::getline(file, line)) {
    auto [expectedCurrentState, expectedSubline] = nextState(state, subline);
    if (state == CSFSParserState::CSFS && expectedCurrentState == CSFSParserState::Time) {
      // moved to new time block, save CSFSentry
      assert(mu > 0);
      assert(from > 0);
      assert(to > from);
      assert(samples > 0);
      parsed.emplace(from, CSFSEntry(timeVector, sizeVector, mu, from, to, samples, csfs));
      mu = from = to = -1.0;
      samples = -1;
    }

    auto nowState = currentState(line);
    if (expectedCurrentState != nowState) {
      throw std::runtime_error(fmt::format("Expected state {}, got {}", expectedCurrentState, currentState));
    }
    std::stringstream tokens(line);
    switch (nowState) {
      case CSFSParserState::Null:
        continue;
      case CSFSParserState::Time:
        tokens >> token;
        while(tokens >> t_dt) timeVector.push_back(t_dt);
        continue;
      case CSFSParserState::Size:
        tokens >> token;
        while(tokens >> t_dt) sizeVector.push_back(t_dt);
        continue;
      case CSFSParserState::Mu:
        tokens >> token; tokens >> mu;
        continue;
      case CSFSParserState::Samples:
        tokens >> token; tokens >> samples;
        continue;
      case CSFSParserState::Interval:
        tokens >> token;
        tokens >> from >> to;
        continue;
      case CSFSParserState::CSFS:
        if (csfs.size() == 0) csfs.resize(3, samples - 1);
        unsigned col = 0;
        while(tokens >> t_dt) csfs(subline, col++) = t_dt;
        continue;
    }
    state = nowState;
    subline = expectedSubline;
  }
  return CSFS(parsed);
}

bool CSFS::verify(std::vector<double> timeVectorOriginal, std::vector<double> sizeVectorOriginal,
    double mu, int samples, std::vector<double> discretizationOriginal) {
  std::vector<double> timeVector(timeVectorOriginal);
  std::vector<double> sizeVector(sizeVectorOriginal);
  std::vector<double> discretization(discretizationOriginal);
  timeVector.pop_back();
  sizeVector.pop_back();
  discretization.pop_back();
  for (double from : discretization) {
    auto search = mCSFS.find(from);
    if (search == mCSFS.end()) {
      fmt::print("Warning:\tCSFS does not contain interval {}.\n", from);
      return false;
    }
    auto thisEntry = mCSFS.at(from);
    if (thisEntry.getMu() != mu) {
      fmt::print("Warning:\tCSFS entry {} has different mu: {}.", from, thisEntry.getMu());
      return false;
    }
    if (thisEntry.getTime() != timeVector) {
      fmt::print("{}\n", thisEntry.getTime());
      fmt::print("{}\n", timeVector);
      fmt::print("Warning:\tCSFS entry {} has different time vector.", from);
      return false;
    }
    if (thisEntry.getSize() != sizeVector) {
      std::cout << "Warning:\tCSFS entry " << from << " has different size vector." << std::endl;
      return false;
    }
    if (thisEntry.getSamples() != samples) {
      // if (samples == Integer.MAX_VALUE)  samples = thisEntry.samples;
      std::cout << "Warning:\tCSFS entry " << from << " has different samples (want: " <<
        samples << ", found: " << thisEntry.getSamples() << ")" << std::endl;
      return false;
    }
  }
  return true;
}

std::string CSFS::toString() const {
  std::string repr;
  for (auto const& x: mCSFS) repr += x.second.toString();
  return repr;
}

void CSFS::fixAscertainment(Data data, int samples, Transition transition) {
    computeArraySamplingFactors(data, samples, transition);
    // CSFS is loaded here, but fixed later.
    for (auto const& entry: mCSFS) mAscertainedCSFS.emplace(entry);
    applyFactors();
    mFoldedAscertainedCSFS = foldCSFS(mAscertainedCSFS);
    mCompressedAscertainedEmissionTable = compressCSFS(mFoldedAscertainedCSFS);
}

mat_dt CSFS::computeClassicEmission(std::vector<double> expectedTimes, double mu) {
  array_dt emissionRow(expectedTimes.size());
  for(unsigned i = 0; i < expectedTimes.size(); i++) {
    emissionRow(i) = std::exp(-2 * expectedTimes[i] * mu);
  }
  mat_dt emission(2, emissionRow.size());
  emission.row(0) = emissionRow;
  emission.row(1) = 1 - emissionRow;
  return emission;
}

void CSFS::computeArraySamplingFactors(Data data, int samples, Transition transition) {
    mSamples = samples;
    auto coalDist = transition.getCoalDist();
    array_dt AFS(samples);
    // double[] AFS = new double[samples];
    // the first entry of the CSFS may not be zero, since it's a shared doubleton
    unsigned counter = 0;
    for (auto &[from, csfsEntry] : mCSFS) {
      auto mat_csfs = csfsEntry.getCSFS();
      for (int row = 0; row < 3; row++) {
        for (int column = 0; column < samples - 1; column++) {
          int pos = row + column;
          if (pos > samples / 2) pos = samples - pos;
          AFS[pos] += coalDist[counter] * mat_csfs(row, column);
        }
      }
      counter++;
    }

    AFS[0] = 0.;
    AFS /= AFS.sum(); // normalize spectrum

    // fold AFS
    int halfTotal = samples / 2;
    for (int i = halfTotal + 1; i < samples; i++) {
        AFS[samples - i] += AFS[i];
        AFS[i] = 0;
    }
    AFS /= AFS.sum(); // normalize spectrum

    // foldedAFS contains probability a site has MAF i given the site is
    // polymorphic in the sequence data
    array_dt foldedAFS = AFS.head(halfTotal + 1);

    // now get foldedAFS_array, the probability a site has MAF i given it is polymorphic in the sample (array)
    mArraySpectrum = ArraySpectrum(data, static_cast<unsigned int>(samples));
    auto foldedAFS_array = mArraySpectrum.getSpectrum();
    mArraySamplingFactors.resize(halfTotal + 1);
    mArraySamplingFactors[0] = 0.0;
    for (unsigned i = 1; i < foldedAFS_array.size(); i++) {
        mArraySamplingFactors[i] = foldedAFS_array[i] / foldedAFS[i];
    }
}

void CSFS::applyFactors() {
  // apply sampling factors and renormalize
  // note that the first entry of the CSFS may not be zero, since it's a shared doubleton
  double monomorphic = mArraySpectrum.getMonomorphic();
  for (auto &[from, csfsEntry] : mAscertainedCSFS) {
    auto thisCSFS = csfsEntry.getCSFS();
    if (thisCSFS.size() > 0) thisCSFS(0, 0) = 0.; else throw std::runtime_error("CSFS is empty!");
    double norm = 0.;
    for (int row = 0; row < 3; row++) {
      for (int column = 0; column < mSamples - 1; column++) {
        // if the spectrum is folded, this emission is mapped to this position
        int pos = row + column;
        if (pos > mSamples / 2) pos = mSamples - pos;
        // and if we're looking at array data, this MAF is adjusted using this factor
        thisCSFS(row, column) *= mArraySamplingFactors[pos];
        // sum value to renomralize to 1 later on
        norm += thisCSFS(row, column);
      }
    }
    norm /= 1 - monomorphic;
    thisCSFS /= norm;
    thisCSFS(0, 0) = monomorphic;
    // ascertainedCSFS.get(from).CSFS = thisCSFS; [notreq, using refs]
    }
}

std::map<double, CSFSEntry> CSFS::foldCSFS(std::map<double, CSFSEntry> csfsMap) {
  std::map<double, CSFSEntry> foldedCSFS;
  int samples = csfsMap.at(0).getSamples();
  int undistinguished = samples - 2;
  for (auto &[from, foldedEntry] : csfsMap) {
    auto thisCsfs_double = foldedEntry.getCSFS();
    // code to fold the spectrum
    if (samples % 2 != 0) throw std::runtime_error("ConditionalSFS called with odd number of samples.");
    int half = samples / 2;
    mat_dt thisCsfs_double_folded(2, half + 1);
    thisCsfs_double_folded.setZero();
    for (int row = 0; row < 3; row++) {
      for (int column = 0; column < undistinguished + 1; column++) {
        auto [dist, undist] = getFoldedObservationFromUnfolded(std::make_pair(row, column), samples);
        thisCsfs_double_folded(dist, undist) += thisCsfs_double(row, column);
      }
    }
    foldedEntry.setCSFS(thisCsfs_double_folded);
    foldedCSFS.emplace(std::make_pair(from, foldedEntry));
    }
  return foldedCSFS;
}

std::pair<int, int> CSFS::getFoldedObservationFromUnfolded(std::pair<int, int> unfolded, int totalSamples) {
  auto [dist, undist] = unfolded;
  if (totalSamples % 2 != 0) throw std::runtime_error(
      "Function getFoldedObservationFromUnfolded was called with odd total sample size. "
      "Only diploid samples are supported at the moment.");
  if (undist + dist > totalSamples / 2) undist = totalSamples - 2 - undist; // flip
  if (dist == 2) dist = 0;
  return std::make_pair(dist, undist);
}

mat_dt CSFS::compressCSFS(std::map<double, CSFSEntry> csfsMap) {
  mat_dt compressed(2, csfsMap.size());
  compressed.setZero();
  int timeInterval = 0;
  for (auto &[from, csfsEntry] : csfsMap) {
    auto thisCSFS = csfsEntry.getCSFS();
    for (int k = 0; k < thisCSFS.cols(); k++) {
      compressed(0, timeInterval) += thisCSFS(0, k);
      compressed(1, timeInterval) += thisCSFS(1, k);
    }
    timeInterval++;
  }
  return compressed;
}

} // namespace asmc
