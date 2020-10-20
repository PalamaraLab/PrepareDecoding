// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#ifndef PREPAREDECODING_TRANSITION_HPP
#define PREPAREDECODING_TRANSITION_HPP

#include <vector>
#include <array>
#include "EigenTypes.hpp"

using asmc::vec_dt;
using asmc::mat_dt;

namespace asmc {

enum TransitionType { SMC, SMC1, CSC };
const std::array<const char*, 3> TransitionType_str = {"SMC", "SMC1", "CSC"};

class Transition {

  private:

    std::vector<double> mTime;
    std::vector<double> mSize;
    std::vector<double> mDiscretization;
    std::vector<double> mTimeVectorPlusInfinity;
    std::vector<double> mExpectedTimes;
    TransitionType mType;
    unsigned int mStates;
    // coalescent arrays that will be computed once and depend on the demography
    std::vector<double> mProbNotCoalesceBetweenExpectedTimes;
    std::vector<double> mProbNotCoalesceBetweenTimeIntervals;
    std::vector<double> mProbCoalesceBetweenExpectedTimesAndUpperLimit;
    vec_dt mColumnRatios;

    static mat_dt identity(TransitionType type);
    static std::tuple<vec_dt, vec_dt, vec_dt, vec_dt> getLinearTimeDecodingQuantitiesGivenTransition(mat_dt T);

    mat_dt transitionMatrix(double r);
    static mat_dt getExponentiatedTransitionMatrix(double N, double r, double time, TransitionType type);
    double getTransitionFromStoInterval(double r, double timeS, double fromTime, double toTime, TransitionType type);
    std::vector<double> expectedIntervalTimesPiecewise();
    double getCumulativeTransitionProbability(double r, double timeS, double timeT, TransitionType type);

    double expectedTimeFromStoT(double timeS, double timeT);
    double coalesceFromStoT(double timeS, double timeT);
    double notCoalesceFromStoT(double timeS, double timeT);
    double getSizeInPiecewiseAtTimeT(double timeT);
    mat_dt computeTransitionPiecewiseUpToTimeT(double r, double time, TransitionType type);
    unsigned int findIntervalForTime(double t);
    mat_dt computeTransitionPiecewiseFromTimeSToTimeT(double r, double timeS, double timeT, TransitionType type);
    double cumulativeCoalesceFromStoTsmart(double timeS, double timeT);
    std::pair<mat_dt, mat_dt> getOmegas(double r, TransitionType type);
    void computeCoalescentVectors();

  public:

    Transition(std::vector<double> timeVector, std::vector<double> sizeVector, std::vector<double> discretization,
               TransitionType type);
    static std::vector<double> getTimeExponentialQuantiles(int numQuantiles, std::vector<double> timeVector,
                                                    std::vector<double> sizeFromVector);
    static std::vector<double> getTimeErlangQuantiles(int numQuantiles, std::vector<double> timeVector,
                                               std::vector<double> sizeFromVector);
    double cumulativeCoalesceFromStoT(double timeS, double timeT);
    std::tuple<vec_dt, vec_dt, vec_dt, vec_dt> getLinearTimeDecodingQuantitiesAndMatrixGivenDistance(double rho);
    std::vector<double>& getTime() { return mTime; }
    std::vector<double>& getSize() { return mSize; }
    std::vector<double>& getDiscretization() { return mDiscretization; }
    std::vector<double>& getExpectedTimes() { return mExpectedTimes; }
    TransitionType getType() { return mType; }
    vec_dt& getColumnRatios() { return mColumnRatios; }

    std::vector<double> getCoalDist();
    static const std::vector<double> EUsize;
    static const std::vector<double> EUtime;
};

} // namespace asmc

#endif // PREPAREDECODING_TRANSITION_HPP
