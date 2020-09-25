// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include <algorithm>
#include <fmt/core.h>
#include "DecodingQuantities.hpp"

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
    if (percentage != lastPercentage) fmt::print("\nGenetic distances progress: {}%", percentage);
    lastPercentage = percentage;
  }
  fmt::print("\n");

  // compute homozygous emissions
  int phys = startPhys;
  while (phys < maxPhys) {
    mPhysDistances.push_back(phys);
    phys = nextPhys(phys);
  }
  for (unsigned i = 0; i < mPhysDistances.size(); i++) {
    mHomozygousEmissions.push_back(CSFS::computeClassicEmission(
          mExpectedTimes, mPhysDistances[i] * mMu));
    int percentage = static_cast<int>(std::round(100 * i / static_cast<double>(mPhysDistances.size())));
    if (percentage != lastPercentage) fmt::print("\rPhysical distances progress: {}%", percentage);
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
  int log10 = static_cast<int>(std::max(0, std::floor(std::log10(phys)) - precision));
  int factor = static_cast<int>(std::pow(10, log10));
  return static_cast<int>(std::round(phys / static_cast<double>(factor + 1)) * factor);
}

double DecodingQuantities::nextGen(double gen) {
  double gene1e10 = gen * 1e10;
  int log10 = static_cast<int>(std::max(0, std::floor(std::log10(gene1e10)) - precision));
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

// save all compute quantities to file
void DecodingQuantities::saveDecodingQuantities(std::string_view outputFileRoot) {
   BufferedWriter outBuff = Utils.openGzipFileForWriting(outputFileRoot + ".decodingQuantities.gz");
   // write model parameters
   outBuff.write("TransitionType\n" + this.transitionType + "\n\n");
   outBuff.write("States\n" + this.states + "\n\n");
   outBuff.write("CSFSSamples\n" + this.CSFSSamples + "\n\n");
   outBuff.write("TimeVector\n");
   outBuff.write(Utils.doubleMatrixToString(timeVector));
   outBuff.write("\n");
   outBuff.write("SizeVector\n");
   outBuff.write(Utils.doubleMatrixToString(sizeVector));
   outBuff.write("\n");
   outBuff.write("Discretization\n");
   outBuff.write(Utils.doubleMatrixToString(discretization));
   outBuff.write("\n");
   outBuff.write("ExpectedTimes\n");
   outBuff.write(Utils.doubleMatrixToString(expectedTimes));
   // write sequence Emissions
   outBuff.write("\n");
   for (int undistinguished = 0; undistinguished < this.CSFSmap.firstEntry().getValue().CSFS[0].length; undistinguished++) {
       outBuff.write("CSFS\t" + undistinguished + "\n");
       for (int distinguished = 0; distinguished < 3; distinguished++) {
           for (double from : this.CSFSmap.keySet()) {
               outBuff.write(this.CSFSmap.get(from).CSFS[distinguished][undistinguished] + "\t");
           }
           outBuff.write("\n");
       }
   }
   outBuff.write("\n");
   for (int undistinguished = 0; undistinguished < this.foldedCSFSmap.firstEntry().getValue().CSFS[0].length; undistinguished++) {
       outBuff.write("FoldedCSFS\t" + undistinguished + "\n");
       for (int distinguished = 0; distinguished < 2; distinguished++) {
           for (double from : this.foldedCSFSmap.keySet()) {
               outBuff.write(this.foldedCSFSmap.get(from).CSFS[distinguished][undistinguished] + "\t");
           }
           outBuff.write("\n");
       }
   }
   outBuff.write("\n");
   outBuff.write("ClassicEmission\n");
   outBuff.write(Utils.doubleMatrixToString(classicEmissionTable));
   // write ascertained Emissions
   outBuff.write("\n");
   for (int undistinguished = 0; undistinguished < this.ascertainedCSFSmap.firstEntry().getValue().CSFS[0].length; undistinguished++) {
       outBuff.write("AscertainedCSFS\t" + undistinguished + "\n");
       for (int distinguished = 0; distinguished < 3; distinguished++) {
           for (double from : this.ascertainedCSFSmap.keySet()) {
               outBuff.write(this.ascertainedCSFSmap.get(from).CSFS[distinguished][undistinguished] + "\t");
           }
           outBuff.write("\n");
       }
   }
   outBuff.write("\n");
   for (int undistinguished = 0; undistinguished < this.foldedAscertainedCSFSap.firstEntry().getValue().CSFS[0].length; undistinguished++) {
       outBuff.write("FoldedAscertainedCSFS\t" + undistinguished + "\n");
       for (int distinguished = 0; distinguished < 2; distinguished++) {
           for (double from : this.foldedAscertainedCSFSap.keySet()) {
               outBuff.write(this.foldedAscertainedCSFSap.get(from).CSFS[distinguished][undistinguished] + "\t");
           }
           outBuff.write("\n");
       }
   }
   outBuff.write("\n");
   outBuff.write("CompressedAscertainedEmission\n");
   outBuff.write(Utils.doubleMatrixToString(compressedEmissionTable));
   // write initial state distribution
   outBuff.write("\n");
   outBuff.write("initialStateProb\n");
   outBuff.write(Utils.doubleMatrixToString(initialStateProb.transpose()));
   // write column ratios
   outBuff.write("\n");
   outBuff.write("ColumnRatios\n");
   outBuff.write(Utils.doubleMatrixToString(columnRatios.transpose()));
   // write row ratios
   outBuff.write("\n");
   outBuff.write("RowRatios\n");
   for (int i = 0; i < rowRatiosVectorMap.size(); i++) {
       DoubleMatrix thisRR = rowRatiosVectorMap.get(i);
       outBuff.write(geneticDistancesList.get(i) + "\t" + Utils.doubleMatrixToString(thisRR.transpose()));
   }
   // write U vectors
   outBuff.write("\n");
   outBuff.write("Uvectors\n");
   for (int i = 0; i < UvectorMap.size(); i++) {
       DoubleMatrix thisU = UvectorMap.get(i);
       outBuff.write(geneticDistancesList.get(i) + "\t" + Utils.doubleMatrixToString(thisU.transpose()));
   }
   // write B vectors
   outBuff.write("\n");
   outBuff.write("Bvectors\n");
   for (int i = 0; i < BvectorMap.size(); i++) {
       DoubleMatrix thisB = BvectorMap.get(i);
       outBuff.write(geneticDistancesList.get(i) + "\t" + Utils.doubleMatrixToString(thisB.transpose()));
   }
   // write D vectors
   outBuff.write("\n");
   outBuff.write("Dvectors\n");
   for (int i = 0; i < DvectorMap.size(); i++) {
       DoubleMatrix thisD = DvectorMap.get(i);
       outBuff.write(geneticDistancesList.get(i) + "\t" + Utils.doubleMatrixToString(thisD.transpose()));
   }
   // write homozygous emissions
   outBuff.write("\n");
   outBuff.write("HomozygousEmissions\n");
   for (int i = 0; i < physDistancesList.size(); i++) {
       int phys = physDistancesList.get(i);
       outBuff.write(phys + "\t" + Utils.doubleMatrixToString(HomozygousEmissionMap.get(i).getRow(0)));
   }
   outBuff.close();
}

} // namespace asmc
