// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include <vector>
#include "EigenTypes.hpp"

#ifndef PREPAREDECODING_TRANSITION_HPP
#define PREPAREDECODING_TRANSITION_HPP

using asmc::vec_dt;
using asmc::mat_dt;

namespace asmc {

enum TransitionType { SMC, SMC1, CSC };

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
    static std::vector<double> getTimeExponentialQuantiles(int numQuantiles, std::vector<double> timeVector,
                                                    std::vector<double> sizeFromVector);
    static std::vector<double> getTimeErlangQuantiles(int numQuantiles, std::vector<double> timeVector,
                                               std::vector<double> sizeFromVector);
    static std::tuple<vec_dt, vec_dt, vec_dt, vec_dt> getLinearTimeDecodingQuantitiesGivenTransition(mat_dt T);

    std::vector<double> getExpectedTimes();
    std::tuple<vec_dt, vec_dt, vec_dt, vec_dt> getLinearTimeDecodingQuantitiesAndMatrixGivenDistance(double rho);
    mat_dt transitionMatrix(double r);
    static mat_dt getExponentiatedTransitionMatrix(double N, double r, double time, TransitionType type);
    double getTransitionFromStoInterval(double r, double timeS, double fromTime, double toTime, TransitionType type);
    std::vector<double> expectedIntervalTimesPiecewise();
    double getCumulativeTransitionPobability(double r, double timeS, double timeT, TransitionType type);

    double expectedTimeFromStoT(double timeS, double timeT);
    double cumulativeCoalesceFromStoT(double timeS, double timeT);
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
    std::vector<double> getCoalDist();
    const std::vector<double> EUsize {
      145041., 129827., 116209., 104020., 93109., 83342., 74600., 66775., 59771., 53501., 47892., 44915., 43684.,
      42486.,  41321.,  40188.,  39086.,  38014., 36972., 35958., 34972., 34013., 33080., 32173., 31291., 30433.,
      29598.,  28787.,  27997.,  27230.,  26483., 25757., 25050., 24364., 23695., 23046., 22414., 21799., 21201.,
      20620.,  20055.,  19505.,  18970.,  18450., 17944., 17452., 16973., 16508., 16055., 15615., 15186., 14770.,
      14365.,  13971.,  13588.,  13215.,  12853., 12501., 12158., 11824., 11500., 11185., 10878., 10580., 10289.,
      10007.,  9733.,   9466.,   9206.,   8954.,  8708.,  8470.,  7871.,  7605.,  6890.,  6080.,  5191.,  4225.,
      3180.,   2843.,   2604.,   2830.,   2940.,  2694.,  5002.,  6138.,  5562.,  6093.,  7825.,  9864.,  13076.,
      15823.,  18475.,  20949.,  23350.,  25032., 26366., 26381., 26035., 24927., 23896., 22998., 23071., 23312.,
      23339.,  23145.,  22510.,  21536.,  20189., 18401., 15099., 15628., 15628.};

    const std::vector<double> EUtime {0., 10., 20., 30., 40., 50., 60., 70., 80.,
      90., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190.,
      200., 210., 220., 230., 240., 250., 260., 270., 280., 290., 300.,
      310., 320., 330., 340., 350., 360., 370., 380., 390., 400., 410.,
      420., 430., 440., 450., 460., 470., 480., 490., 500., 510., 520.,
      530., 540., 550., 560., 570., 580., 590., 600., 610., 620., 630.,
      640., 650., 660., 670., 680., 690., 700., 710., 720., 724., 795.,
      874., 960., 1055., 1160., 1275., 1402., 1541., 1694., 1862., 2047.,
      2250., 2474., 2720., 2990., 3287., 3614., 3973., 4368., 4802., 5280.,
      5805., 6382., 7017., 7715., 8482., 9325., 10252., 11271., 12391.,
      13623., 14977., 16466., 18102., 19901., 21879., 24053., 26443.,
      std::numeric_limits<double>::infinity()};
};

#endif // PREPAREDECODING_TRANSITION_HPP

} // namespace asmc
