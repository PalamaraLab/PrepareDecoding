// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include <algorithm>
#include <fstream>
#include <iostream>
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <fmt/ostream.h>
#include "DecodingQuantities.hpp"
#include "Utils.hpp"

namespace asmc {

DecodingQuantities::DecodingQuantities(CSFS& csfs, Transition& transition, double mu) {
  mTime = transition.getTime();
  mMu = mu;
  mSize = transition.getSize();
  mCSFSSamples = csfs.getSamples();
  mDiscretization = transition.getDiscretization();
  mExpectedTimes = transition.getExpectedTimes();
  mStates = static_cast<unsigned int>(mDiscretization.size()) - 1;
  mColumnRatios = transition.getColumnRatios();
  mTransitionType = transition.getType();
  mCSFS = csfs.getCSFS();
  mFoldedCSFS = csfs.getFoldedCSFS();
  mClassicEmissionTable = CSFS::computeClassicEmission(mExpectedTimes, mMu);
  mAscertainedCSFS = csfs.getAscertainedCSFS();
  mFoldedAscertainedCSFS = csfs.getFoldedAscertainedCSFS();
  mCompressedEmissionTable = csfs.getCompressedAscertainedEmissionTable();
  int lastPercentage = -1;

  // compute transition quantities
  // putting 0 in. This is not really used, depending on the value of smallNumber in the decoding side.
  computeTransitionQuantitiesForOnePosition(0.0, transition);
  mGeneticDistances.push_back(0.0);
  double genDist = startGen;
  while (genDist < maxGen) {
    mGeneticDistances.push_back(genDist);
    genDist = nextGen(genDist);
  }
        // start at 1, 0 is done already
  for (unsigned i = 1; i < mGeneticDistances.size(); i++) {
    genDist = mGeneticDistances[i];
    // compute emission and transition quantities in linear time
    computeTransitionQuantitiesForOnePosition(genDist, transition);
    int percentage = static_cast<int>(
        std::round(100 * i / static_cast<double>(mGeneticDistances.size())));
    if (percentage != lastPercentage)
      std::cout << "Genetic distances progress: " << percentage << "%\t\r" << std::flush;
    lastPercentage = percentage;
  }
  fmt::print("\n");

  // compute homozygous emissions
  int phys = startPhys;
  while (phys < maxPhys) mPhysDistances.push_back(phys = nextPhys(phys));
  for (unsigned i = 0; i < mPhysDistances.size(); i++) {
    mHomozygousEmissions.push_back(CSFS::computeClassicEmission(
          mExpectedTimes, mPhysDistances[i] * mMu));
    int percentage = static_cast<int>(std::round(100 * i / static_cast<double>(mPhysDistances.size())));
    if (percentage != lastPercentage)
      std::cout << "Physical distances progress: " << percentage << "%\t\r" << std::flush;
    lastPercentage = percentage;
  }
  fmt::print("\n");
   // get initial state distribution from coalescent distribution
  double lastProb = 0;
  vec_dt cumProb(mStates);
  mInitialStateProb.resize(mStates);
  for (unsigned i = 0; i < mStates; i++) {
    double timeT = mDiscretization.at(i + 1);
    cumProb[i] = transition.cumulativeCoalesceFromStoT(0, timeT);
    mInitialStateProb[i] = cumProb[i] - lastProb;
    lastProb = cumProb[i];
  }

}

int DecodingQuantities::nextPhys(int phys) {
  if (phys < 0) throw std::runtime_error("Int overflow: " + std::to_string(phys));
  int log10 = static_cast<int>(std::max(0.0, std::floor(std::log10(phys)) - precision));
  int factor = static_cast<int>(std::pow(10, log10));
  return static_cast<int>(std::round(phys / static_cast<double>(factor) + 1) * factor);
}

double DecodingQuantities::nextGen(double gen) {
  double gene1e10 = gen * 1e10;
  int log10 = static_cast<int>(std::max(0.0, std::floor(std::log10(gene1e10)) - precision));
  double factor = std::pow(10, log10);
  return (std::round(gene1e10 / factor) + 1) * factor / 1e10;
}

// compute quantities for a give genetic distance
void DecodingQuantities::computeTransitionQuantitiesForOnePosition(double genDist, Transition& transition) {
    auto [D, B, U, RR] = transition.getLinearTimeDecodingQuantitiesAndMatrixGivenDistance(genDist);
    mDvectors.push_back(D);
    mBvectors.push_back(B);
    mUvectors.push_back(U);
    mRowRatioVectors.push_back(RR);
}

std::string DecodingQuantities::csfsToString(const std::string& header, const std::map<double, CSFSEntry>& csfs_map) {
  std::string out;
  for (int undistinguished = 0; undistinguished < (csfs_map.begin()->second).getCSFS().cols(); undistinguished++) {
    out += fmt::format("{}\t{}\n", header, undistinguished);
    for (int distinguished = 0; distinguished < 3; distinguished++) {
      for (auto const& [from, csfs] : csfs_map) {
        out += fmt::format("{}\t", csfs.getCSFS()(distinguished, undistinguished));
      }
      out += "\n";
    }
   }
  return out + "\n";
}

std::string DecodingQuantities::vectorsToString(const std::string& header, const std::vector<vec_dt>& vector_list) {
  std::string out(header + "\n");
  for (unsigned i = 0; i < vector_list.size(); i++) {
    out += fmt::format("{}\t{}\n", mGeneticDistances.at(i), vector_list.at(i));
  }
  return out + "\n";
}

// save all compute quantities to file
void DecodingQuantities::saveDecodingQuantities(std::string_view outputFileRoot) {
  gzFile file = gzopen(fmt::format("{}.decodingQuantities.gz", outputFileRoot).c_str(), "wb");
  // write model parameters
  writegz(file, fmt::format("TransitionType\n{}\n\n", mTransitionType));
  writegz(file, fmt::format("States\n{}\n\n", mStates));
  writegz(file, fmt::format("CSFSSamples\n{}\n\n", mCSFSSamples));
  writegz(file, fmt::format("TimeVector\n{}\nSizeVector\n{}\nDiscretization\n{}\nExpectedTimes\n{}\n",
        mTime, mSize, mDiscretization, mExpectedTimes));
  // write sequence Emissions
  writegz(file, csfsToString("CSFS", mCSFS));
  writegz(file, csfsToString("FoldedCSFS", mFoldedCSFS));
  writegz(file, fmt::format("ClassicEmission\n{}\n", mClassicEmissionTable));
  writegz(file, csfsToString("AscertainedCSFS", mAscertainedCSFS));
  writegz(file, csfsToString("FoldedAscertainedCSFS", mFoldedAscertainedCSFS));
  writegz(file, fmt::format("CompressedAscertainedEmission\n{}\n", mCompressedEmissionTable));
   // write initial state distribution
  writegz(file, fmt::format("initialStateProb\n{}\n", mInitialStateProb));
  writegz(file, fmt::format("ColumnRatios\n{}\n", mColumnRatios));
  writegz(file, vectorsToString("RowRatios", mRowRatioVectors));
  writegz(file, vectorsToString("Uvectors", mUvectors));
  writegz(file, vectorsToString("Bvectors", mBvectors));
  writegz(file, vectorsToString("Dvectors", mDvectors));
  // write homozygous emissions
  writegz(file, "HomozygousEmissions\n");
  for (unsigned i = 0; i < mPhysDistances.size(); i++) {
    auto phys = mPhysDistances.at(i);
    writegz(file, fmt::format("{}\t{}\n", phys, mHomozygousEmissions.at(i).row(0)));
   }
   gzclose(file);
}

void DecodingQuantities::saveIntervals(std::string_view outputFileRoot) {
  std::ofstream file(fmt::format("{}.intervalsInfo", outputFileRoot));
  for (unsigned i = 0; i < mExpectedTimes.size(); i++)
    file << mDiscretization.at(i) << "\t" << mExpectedTimes.at(i) <<
    "\t" << mDiscretization.at(i + 1) << std::endl;
  file.close();
}

} // namespace asmc
