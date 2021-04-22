// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#ifndef PREPAREDECODING_DECODINGQUANTITIES_HPP
#define PREPAREDECODING_DECODINGQUANTITIES_HPP

#include <vector>
#include "EigenTypes.hpp"
#include "Csfs.hpp"
#include "Transition.hpp"

using asmc::vec_dt;
using asmc::mat_dt;

namespace asmc {

class DecodingQuantities {

  private:
    // these are the quantities loaded from file, which are required for the linear-time decoding
    TransitionType mTransitionType;
    vec_dt mInitialStateProb;
    std::vector<double> mTime;
    std::vector<double> mSize;
    std::vector<double> mDiscretization;
    std::vector<double> mExpectedTimes;
    std::vector<double> mGeneticDistances;
    std::vector<int> mPhysDistances;
    std::vector<vec_dt> mDvectors;
    std::vector<vec_dt> mBvectors;
    std::vector<vec_dt> mUvectors;
    std::vector<vec_dt> mRowRatioVectors;
    std::vector<mat_dt> mHomozygousEmissions;
    vec_dt mColumnRatios;

    unsigned int mStates, mCSFSSamples;
    double mMu;

    // sequence emission
    std::map<double, CSFSEntry> mCSFS;
    std::map<double, CSFSEntry> mFoldedCSFS;
    mat_dt mClassicEmissionTable;

    // ascertained emission
    std::map<double, CSFSEntry> mAscertainedCSFS;
    std::map<double, CSFSEntry> mFoldedAscertainedCSFS;
    mat_dt mCompressedEmissionTable;

    // for a number between 10^(n) and 10^(n+1), proceed in steps of size max(1, 10^(n-precision))
    constexpr static int precision = 2;
    constexpr static double minGenetic = 1e-10;

    // start/end for physical and genetic distances (in Morgans)
    constexpr static double startGen = 1e-10;
    constexpr static double maxGen = 0.3; // 30 centiMorgans
    constexpr static int startPhys = 1;
    constexpr static int maxPhys = 100000000; // 100 Mb

    void computeTransitionQuantitiesForOnePosition(double genDist, Transition& transition);
    static std::string csfsToString(const std::string& header, const std::map<double, CSFSEntry>& csfs_map);
    std::string vectorsToString(const std::string& header, const std::vector<vec_dt>& vector_list);

  public:

    // for each genetic distance, compute decoding quantities
    DecodingQuantities(CSFS& csfs, Transition& transition, double mu);
    static int nextPhys(int phys);
    static double nextGen(double gen);

    /**
     * Save decoding quantities to file
     * @param outputFileRoot the decoding quantities file will be saved as {outputFileRoot}.decodingQuantities.gz
     */
    void saveDecodingQuantities(std::string_view outputFileRoot);

    /**
     * Save intervals to file
     * @param outputFileRoot the intervals file will be saved as {outputFileRoot}.intervalsInfo
     */
    void saveIntervals(std::string_view outputFileRoot);

    /**
     * Save CSFS to file
     * @param outputFileRoot the CSFS file will be saved as {outputFileRoot}.csfs
     */
    void saveCsfs(std::string_view outputFileRoot);

    TransitionType getTransitionType() { return mTransitionType; }
    vec_dt getInitialStateProb() const { return mInitialStateProb; }
    std::vector<double> getTimes() const { return mTime; }
    std::vector<double> getSizes() const { return mSize; }
    std::vector<double> getDiscretization() const { return mDiscretization; }
    std::vector<double> getExpectedTimes() const { return mExpectedTimes; }
    std::vector<double> getGeneticDistances() const { return mGeneticDistances; }
    std::vector<int> getPhysDistances() const { return mPhysDistances; }
    std::vector<vec_dt> getDvectors() const { return mDvectors; }
    std::vector<vec_dt> getBvectors() const { return mBvectors; }
    std::vector<vec_dt> getUvectors() const { return mUvectors; }
    std::vector<vec_dt> getRowRatioVectors() const { return mRowRatioVectors; }
    std::vector<mat_dt> getHomozygousEmissions() const { return mHomozygousEmissions; }
    vec_dt getColumnRatios() const { return mColumnRatios; }
    unsigned int getStates() const { return mStates; }
    unsigned int getCSFSSamples() const { return mCSFSSamples; }
    double getMu() const { return mMu; }
    std::map<double, CSFSEntry> getCSFS() { return mCSFS; }
    std::map<double, CSFSEntry> getFoldedCSFS() { return mFoldedCSFS; }
    std::map<double, CSFSEntry> getAscertainedCSFS() { return mAscertainedCSFS; }
    std::map<double, CSFSEntry> getFoldedAscertainedCSFS() { return mFoldedAscertainedCSFS; }
    mat_dt getClassicEmission() const { return mClassicEmissionTable; }
    mat_dt getCompressedEmission() const { return mCompressedEmissionTable; }

};

} // namespace asmc

#endif // PREPAREDECODING_DECODINGQUANTITIES_HPP
