// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include "Csfs.hpp"

#include "EigenTypes.hpp"

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <stdexcept>

namespace asmc {


CSFS::CSFS(std::map<double, CSFSEntry> CSFS_) : mCSFS(CSFS_),
  mArraySpectrum({}), mCSFS({}), mAscertainedCSFS({}),
  mFoldedAscertainedCSFS({}), mCompressedAscertainedEmissionTable({}),
  mArraySamplingFactors({}), mSamples(0) {
  if (mCSFS.size() > 0) mFoldedCSFS = foldCSFS(mCSFS);
}

CSFSParserState CSFS::currentState(const std::string& line) {
  std::stringstream lineStream(line);
  std::string token;
  lineStream >> token;
  auto search = stateMap.find(token);
  return ((search != stateMap.end()) ? search->second : CSFSParserState::CSFS);
}

CSFSParserState CSFS::nextState(CSFSParserState state) {
  switch (state) {
    case CSFSParserState::Null:
      return CSFSParserState::Time;

    case CSFSParserState::Time;
      return CSFSParserState::Size;

    case CSFSParserState::Size;
      return CSFSParserState::Mu;

    case CSFSParserState::Mu;
      return CSFSParserState::Samples;

    case CSFSParserState::Samples;
      return CSFSParserState::Interval;

    case CSFSParserState::Interval;
      return CSFSParserState::CSFS;

    case CSFSParserState::CSFS;
      if (mCSFSLines == 3) return CSFSParserState::Time;
      mCSFSLines++;
      return CSFSParserState::CSFS;
  }
}

CSFS CSFS::loadFromFile(std::string_view filename) {
  std::ifstream file(filename);
  std::string line;
  std::map<double, CSFSEntry> parsed;
  std::vector<double> timeVector = {};
  std::vector<double> sizeVector = {},
  mat_dt csfs;
  std::string token;
  double t_dt = {};
  CSFSParserState state = CSFSParserState::Null;
  while(std::getline(file, line)) {
    auto expectedCurrentState = nextState(state);
    if (state == CSFSParserState::CSFS && expectedCurrentState == CSFSParserState::Time) {
      // moved to new time block, save CSFSentry
      assert(mu > 0);
      assert(from > 0);
      assert(to > from);
      assert(samples > 0);
      parsed.emplace(from, CSFSEntry(timeVector, sizeVector, mu, from, to, samples, CSFS));
      mu = from = to = -1.0;
      samples = -1;
    }

    auto currentState = getCurrentState(line);
    if (expectedCurrentState != currentState) {
      throw std::runtime_error(fmt::format("Expected state {}, got {}", expectedCurrentState, currentState));
    }
    std::stringstream tokens(line);
    switch (currentState) {
      case CSFSParserState::Time:
        tokens >> token;
        while(tokens >> t_dt) timeVector.push_back(t_dt);
        continue;
      case CSFSParserState::Size:
        tokens >> token;
        while(tokens >> t_dt) sizeVector.push_back(t_dt);
        continue;
      case CSFSParserState::Mu:
        tokens >> token;
        double mu; tokens >> mu;
        continue;
      case CSFSParserState::Samples:
        tokens >> token;
        int samples; tokens >> samples;
        continue;
      case CSFSParserState::Interval:
        tokens >> token;
        double from, to; tokens >> from >> to;
        continue;
      case CSFSParserState::CSFS:
        if (csfs.size() == 0) csfs.resize(3, samples - 1);
        unsigned col = 0;
        while(tokens >> t_dt) csfs(mCSFSLine, col++) = t_dt;
        continue;
    }
    state = currentState;
  }
  mCSFSLine = 0;
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
      std::cout << "Warning:\tCSFS does not contain interval " << from << "." << std::endl;
      return false;
    }
    auto thisEntry = mCSFS[from];
    if (thisEntry.getMu() != mu) {
      std::cout << "Warning:\tCSFS entry " << from << " has different mu: " << thisEntry.getMu() << "." << std::endl;
      return false;
    }
    if (thisEntry.getTime() != timeVector) {
      std::cout << thisEntry.getTime() << std::endl;
      std::cout << timeVector << std::endl;
      std::cout << "Warning:\tCSFS entry " << from << " has different time vector." << std::endl;
      return false;
    }
    if (thisEntry.getSize() != sizeVector) {
      std::cout << "Warning:\tCSFS entry " << from << " has different size vector." << std::endl;
      return false;
    }
    if (thisEntry.getSamples() != samples) {
      // if (samples == Integer.MAX_VALUE)  samples = thisEntry.samples;
      std::cout << "Warning:\tCSFS entry " << from << " has different samples (want: " <<
        samples << ", found: " << thisEntry.samples << ")" << std::endl;
      return false;
    }
  }
  return true;
}

std::string CSFS::toString() {
  std::string repr;
  for (auto const& x: mCSFS) repr += x.second.toString();
  return repr;
}

void CSFS::fixAscertainment(Data data, int samples, Transition transition) {
    computeArraySamplingFactors(data, samples, transition);
    // CSFS is loaded here, but fixed later.
    for (std::pair<double, CSFSEntry> entry: mCSFS) mAscertainedCSFS.emplace(entry);
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
    std::vector<double> coalDist = transition.getCoalDist();
    std::vector<double> AFS;
    // double[] AFS = new double[samples];
    // the first entry of the CSFS may not be zero, since it's a shared doubleton
    int counter = 0;
    for (auto const& [from, csfs] : mCSFS) {
        double coalescentProbabilityThisTimeSlice = coalDist[counter];
        CSFSEntry thisCsfsDM = CSFS.get(from);
        for (int row = 0; row < 3; row++) {
            for (int column = 0; column < samples - 1; column++) {
                int pos = row + column;
                if (pos > samples / 2) {
                    pos = samples - pos;
                }
                AFS[pos] += coalescentProbabilityThisTimeSlice * thisCsfsDM.CSFS[row][column];
            }
        }
        counter++;
    }

    // normalize spectrum
    AFS[0] = 0.;
    double norm = 0.;
    for (int i = 0; i < samples; i++) {
        norm += AFS[i];
    }
    for (int i = 0; i < samples; i++) {
        AFS[i] /= norm;
    }

    // fold AFS
    int halfTotal = samples / 2;
    for (int i = halfTotal + 1; i < samples; i++) {
        AFS[samples - i] += AFS[i];
        AFS[i] = 0;
    }
    // normalize spectrum
    norm = 0.;
    for (int i = 0; i < samples; i++) {
        norm += AFS[i];
    }
    for (int i = 0; i < samples; i++) {
        AFS[i] /= norm;
    }

    // foldedAFS contains probability a site has MAF i given the site is polymorphic in the sequence data
    double[] foldedAFS = new double[halfTotal + 1];
    for (int i = 0; i <= halfTotal; i++) {
        foldedAFS[i] = AFS[i];
    }

    // now get foldedAFS_array, the probability a site has MAF i given it is polymorphic in the sample (array)
    this.arraySpectrum = new ArraySpectrum(data, samples);
    double[] foldedAFS_array = arraySpectrum.spectrum;
    for (int i = 0; i < foldedAFS_array.length; i++) {
    }
    double[] samplingFactors = new double[halfTotal + 1];

    for (int i = 1; i < foldedAFS_array.length; i++) {
        samplingFactors[i] = foldedAFS_array[i] / foldedAFS[i];
    }
    arraySamplingFactors = samplingFactors;
}

private void applyFactors() {
    // apply sampling factors and renormalize
    // note that the first entry of the CSFS may not be zero, since it's a shared doubleton
    for (double from : ascertainedCSFS.keySet()) {
        double[][] thisCSFS = ascertainedCSFS.get(from).CSFS;
        thisCSFS[0][0] = 0.;
        double norm = 0.;
        for (int row = 0; row < 3; row++) {
            for (int column = 0; column < samples - 1; column++) {
                // if the spectrum is folded, this emission is mapped to this position
                int pos = row + column;
                if (pos > samples / 2) {
                    pos = samples - pos;
                }
                // and if we're looking at array data, this MAF is adjusted using this factor
                double factor = arraySamplingFactors[pos];
                // apply factor
                thisCSFS[row][column] *= factor;
                // sum value to renomralize to 1 later on
                norm += thisCSFS[row][column];
            }
        }
        norm /= 1 - arraySpectrum.monomorphic;
        for (int row = 0; row < 3; row++) {
            for (int column = 0; column < samples - 1; column++) {
                thisCSFS[row][column] /= norm;
            }
        }
        thisCSFS[0][0] = arraySpectrum.monomorphic;
        ascertainedCSFS.get(from).CSFS = thisCSFS;
    }
}

public TreeMap<Double, CSFSEntry> foldCSFS(TreeMap<Double, CSFSEntry> CSFS) {
    TreeMap<Double, CSFSEntry> foldedCSFS = new TreeMap<Double, CSFSEntry>();
    int samples = CSFS.firstEntry().getValue().samples;
    int undistinguished = samples - 2;
    for (double from : CSFS.keySet()) {
        CSFSEntry foldedEntry = new CSFSEntry(CSFS.get(from));
        double[][] thisCsfs_double = foldedEntry.CSFS;
        // code to fold the spectrum
        if (samples % 2 != 0) {
            Utils.exit("ConditionalSFS called with odd number of samples.");
        }
        int half = samples / 2;
        double[][] thisCsfs_double_folded = new double[2][half + 1];
        for (int row = 0; row < 3; row++) {
            for (int column = 0; column < undistinguished + 1; column++) {
                Integer[] coord = new Integer[]{row, column};
                Integer[] foldedCoord = getFoldedObservationFromUnfolded(coord, samples);
                thisCsfs_double_folded[foldedCoord[0]][foldedCoord[1]] += thisCsfs_double[row][column];
            }
        }
        foldedEntry.CSFS = thisCsfs_double_folded;
        foldedCSFS.put(from, foldedEntry);
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

public double[][] compressCSFS(TreeMap<Double, CSFSEntry> CSFS) {
    double[][] compressed = new double[2][CSFS.size()];
    int timeInterval = 0;
    for (double from : CSFS.keySet()) {
        double[][] thisCSFS = CSFS.get(from).CSFS;
        for (int k = 0; k < thisCSFS[0].length; k++) {
            compressed[0][timeInterval] += thisCSFS[0][k];
            compressed[1][timeInterval] += thisCSFS[1][k];
        }
        timeInterval++;
    }
    return compressed;
}

} // namespace asmc
