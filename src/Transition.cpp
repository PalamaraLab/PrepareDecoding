//    This file is part of ASMC, developed by Pier Francesco Palamara.
//
//    ASMC is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    ASMC is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with ASMC.  If not, see <https://www.gnu.org/licenses/>.

#include <vector>
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <unsupported/Eigen/MatrixFunctions>
#include "EigenTypes.hpp"
#include "Transition.hpp"

namespace asmc {
    
const std::vector<double> Transition::EUsize{ 145041., 129827., 116209.,
  104020., 93109., 83342., 74600., 66775., 59771., 53501., 47892., 44915.,
  43684., 42486.,  41321.,  40188.,  39086.,  38014., 36972., 35958.,
  34972., 34013., 33080., 32173., 31291., 30433., 29598.,  28787.,  27997.,
  27230.,  26483., 25757., 25050., 24364., 23695., 23046., 22414., 21799.,
  21201., 20620.,  20055.,  19505.,  18970.,  18450., 17944., 17452.,
  16973., 16508., 16055., 15615., 15186., 14770., 14365.,  13971.,  13588.,
  13215.,  12853., 12501., 12158., 11824., 11500., 11185., 10878., 10580.,
  10289., 10007.,  9733.,   9466.,   9206.,   8954.,  8708.,  8470.,
  7871.,  7605.,  6890.,  6080.,  5191.,  4225., 3180.,   2843.,   2604.,
  2830.,   2940.,  2694.,  5002.,  6138.,  5562.,  6093.,  7825.,  9864.,
  13076., 15823.,  18475.,  20949.,  23350.,  25032., 26366., 26381.,
  26035., 24927., 23896., 22998., 23071., 23312., 23339.,  23145.,  22510.,
  21536.,  20189., 18401., 15099., 15628., 15628.};

const std::vector<double> Transition::EUtime{0., 10., 20., 30., 40., 50.,
  60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180.,
  190., 200., 210., 220., 230., 240., 250., 260., 270., 280., 290., 300.,
  310., 320., 330., 340., 350., 360., 370., 380., 390., 400., 410., 420.,
  430., 440., 450., 460., 470., 480., 490., 500., 510., 520., 530., 540.,
  550., 560., 570., 580., 590., 600., 610., 620., 630., 640., 650., 660.,
  670., 680., 690., 700., 710., 720., 724., 795., 874., 960., 1055., 1160.,
  1275., 1402., 1541., 1694., 1862., 2047., 2250., 2474., 2720., 2990.,
  3287., 3614., 3973., 4368., 4802., 5280., 5805., 6382., 7017., 7715.,
  8482., 9325., 10252., 11271., 12391., 13623., 14977., 16466., 18102.,
  19901., 21879., 24053., 26443.};

std::vector<double> Transition::getTimeExponentialQuantiles(
  int numQuantiles, std::vector<double> timeVector,
  std::vector<double> sizeFromVector) {

  double slice = 1. / numQuantiles;
  double nextQuant = slice;
  double timeStep = 0.1;
  std::vector<double> quantiles{0.};
  double pNotCoal = 1.;
  for (size_t i = 0; i < timeVector.size() - 1; i++) {
    double notCoalRate = 1 - timeStep / sizeFromVector[i];
    for (double t = timeVector[i]; t < timeVector[i + 1]; t += timeStep) {
      pNotCoal *= notCoalRate;
      if (1 - pNotCoal > nextQuant) {
        nextQuant += slice;
        quantiles.push_back(std::round(t * 1000.) / 1000.);
        if (nextQuant >= 1.0 - 1E-10) return quantiles;
      }
    }
  }
  return quantiles;
}

    // unoptimized
std::vector<double> Transition::getTimeErlangQuantiles(int numQuantiles, std::vector<double> timeVector,
                                           std::vector<double> sizeFromVector) {
  double slice = 1. / numQuantiles;
  double nextQuant = slice;
  double timeStep = 0.1;
  std::vector<double> quantiles;
  quantiles.push_back(0.);
  double normalizer = 0.;
  double pNotCoal = 1.;
  const double MAX_T = sizeFromVector.back() * 20; //  10 times the ancestral size
  for (unsigned i = 0; i < timeVector.size() - 1; i++) {
    double coalRate = timeStep / sizeFromVector[i];
    double notCoalRate = 1 - coalRate;
    for (double t = timeVector[i]; t < timeVector[i + 1] && t < MAX_T; t += timeStep) {
      pNotCoal *= notCoalRate;
      normalizer += t * coalRate * pNotCoal;
    }
  }
  pNotCoal = 1.;
  double normalizedCumulative = 0.;
  for (unsigned i = 0; i < timeVector.size() - 1; i++) {
    double coalRate = timeStep / sizeFromVector[i];
    double notCoalRate = 1 - coalRate;
    for (double t = timeVector[i]; t < timeVector[i + 1] && t < MAX_T; t += timeStep) {
      pNotCoal *= notCoalRate;
      normalizedCumulative += t * coalRate * pNotCoal / normalizer;
      if (normalizedCumulative >= nextQuant) {
        nextQuant += slice;
        quantiles.push_back(std::round(t * 1000.) / 1000.);
        if (nextQuant >= 1.0) return quantiles;
      }   
    }
  }
  return quantiles;
}

Transition::Transition(std::vector<double> timeVector, std::vector<double> sizeVector, std::vector<double> discretization, TransitionType type) :
  mTime(std::move(timeVector)), mSize(std::move(sizeVector)), mDiscretization(std::move(discretization)),
  mType(type) {
  // Logger.getLogger().setLevel(LOGLEVEL)i;
  mTimeVectorPlusInfinity = mTime;
  mTimeVectorPlusInfinity.push_back(std::numeric_limits<double>::infinity());
  mExpectedTimes = expectedIntervalTimesPiecewise();
  mStates = static_cast<unsigned int>(mDiscretization.size()) - 1;
  computeCoalescentVectors();
}

mat_dt Transition::identity(TransitionType type) {
  if (type == CSC) {
    return four_dt::Identity();
  } else {
    return three_dt::Identity();
  }
}
    
std::tuple<vec_dt, vec_dt, vec_dt, vec_dt> Transition::getLinearTimeDecodingQuantitiesAndMatrixGivenDistance(double rho) {
  auto [omegasAtBoundaries, omegasAtExpectedTimes] = getOmegas(rho, mType);
  vec_dt D(mStates);
  vec_dt B(mStates - 1);
  vec_dt U(mStates - 1);
  vec_dt RR(mStates - 1);
  double diagonal;
  vec_dt omegaAtBoundaries, omegaAtExpectedTimes;
    /* DoubleMatrix D = new DoubleMatrix(states); */
    /* DoubleMatrix B = new DoubleMatrix(states - 1); */
    /* DoubleMatrix U = new DoubleMatrix(states - 1); */
    /* DoubleMatrix RR = new DoubleMatrix(states - 1); */

  // others should be computed only for i > 0
  omegaAtBoundaries = omegasAtBoundaries.row(0);
  omegaAtExpectedTimes = omegasAtExpectedTimes.row(0);
  diagonal = (mType == CSC)
    ? omegaAtExpectedTimes(0, 0) + mProbCoalesceBetweenExpectedTimesAndUpperLimit[0] *
      (omegaAtExpectedTimes(0, 1) + omegaAtExpectedTimes(0, 2)) + omegaAtExpectedTimes(0, 3) - omegaAtBoundaries(0, 3)
    : omegaAtExpectedTimes(0, 0) + mProbCoalesceBetweenExpectedTimesAndUpperLimit[0] *
      omegaAtExpectedTimes(0, 1) + omegaAtExpectedTimes(0, 2) - omegaAtBoundaries(0, 2);
  D[0] = diagonal;

  // now compute all for each i
  for (unsigned i = 1; i < mStates; i++) {
    omegaAtBoundaries = omegasAtBoundaries.row(i);
    omegaAtExpectedTimes = omegasAtExpectedTimes.row(i);
    diagonal = (mType == CSC)
      ? omegaAtExpectedTimes(0, 0) + mProbCoalesceBetweenExpectedTimesAndUpperLimit[i] *
        (omegaAtExpectedTimes(0, 1) + omegaAtExpectedTimes(0, 2)) + omegaAtExpectedTimes(0, 3) - omegaAtBoundaries(0, 3)
      : omegaAtExpectedTimes(0, 0) + mProbCoalesceBetweenExpectedTimesAndUpperLimit[i] *
        omegaAtExpectedTimes(0, 1) + omegaAtExpectedTimes(0, 2) - omegaAtBoundaries(0, 2);
    D[i] = diagonal;
    B[i - 1] = ((mType == CSC)
      ? omegaAtBoundaries(0, 3) - omegasAtBoundaries.row(i - 1)(0, 3)
      : omegaAtBoundaries(0, 2) - omegasAtBoundaries.row(i - 1)(0, 2));
  }
  // do U and RR up to states - 2
  for (unsigned i = 0; i < mStates - 2; i++) {
    double omegaSi = (mType == CSC)
      ? omegasAtExpectedTimes.row(i)(0, 1) + omegasAtExpectedTimes.row(i)(0, 2)
      : omegasAtExpectedTimes.row(i)(0, 1);
    double omegaSiplus1 = (mType == CSC)
      ? omegasAtExpectedTimes.row(i + 1)(0, 1) + omegasAtExpectedTimes.row(i + 1)(0, 2)
      : omegasAtExpectedTimes.row(i + 1)(0, 1);
    U[i] = omegaSi * (1 - mProbCoalesceBetweenExpectedTimesAndUpperLimit[i]) *
      (1 - mProbNotCoalesceBetweenTimeIntervals[i + 1]);
    // rho == 0 --> transition is identity, ratios are 0/0. Avoid NaN by setting to 1.
    RR[i] = ((rho == 0) ? 1. : omegaSi * mProbNotCoalesceBetweenExpectedTimes[i] / omegaSiplus1);
  }
  // do last for U
  double omegaSi = (mType == CSC)
    ? omegasAtExpectedTimes.row(mStates - 2)(0, 1) + omegasAtExpectedTimes.row(mStates - 2)(0, 2)
    : omegasAtExpectedTimes.row(mStates - 2)(0, 1);
  U[mStates - 2] = omegaSi * (1 - mProbCoalesceBetweenExpectedTimesAndUpperLimit[mStates - 2]) *
    (1 - mProbNotCoalesceBetweenTimeIntervals[mStates - 1]);

  // TODO: What goes in last element of RR
  return std::make_tuple(D, B, U, RR);
}

// note: can compute these in linear time instead of building transition matrix (quadratic)
std::tuple<vec_dt, vec_dt, vec_dt, vec_dt> Transition::getLinearTimeDecodingQuantitiesGivenTransition(mat_dt T) {
  vec_dt D;
  auto N = static_cast<unsigned int>(T.cols());
  vec_dt B(N - 1);
  vec_dt U(N - 1);
  vec_dt RR(N - 1);
  D = T.diagonal();
  for (unsigned i = 0; i < N - 1; i++) B[i] = T(i + 1, i); // below
  for (unsigned i = 0; i < N - 1; i++) U[i] = T(i, i + 1); // above
  for (unsigned i = 0; i < N - 2; i++) {
    // ratio of columns
    if (T(i, N - 1) == T(i + 1, N - 1))
      RR[i] = 1.; // avoids 0/0 for totally linked sites in which T = identity
    else
      RR[i] = T(i, N - 1) / T(i + 1, N - 1);
  }
  return std::make_tuple(D, B, U, RR);
}

mat_dt Transition::transitionMatrix(double r) {
  std::vector<double> expIntervals = expectedIntervalTimesPiecewise();
  mat_dt matrix;
  int nExpIntervals = static_cast<int>(expIntervals.size());
  matrix.resize(nExpIntervals, nExpIntervals);
  for (unsigned i = 0; i < mDiscretization.size() - 1; i++) {
    double timeS = expIntervals[i];
    for (unsigned j = 0; j < mDiscretization.size() - 1; j++) {
      double fromTime = mDiscretization[j];
      double toTime = mDiscretization[j + 1];
      matrix(i, j) = getTransitionFromStoInterval(r, timeS, fromTime, toTime, mType);
    }
  }
  return matrix;
}

// can change this for other models (e.g. piecewise exponential)
mat_dt Transition::getExponentiatedTransitionMatrix(double N, double r, double time, TransitionType type) {
  three_dt mat3;
  double rho = 2 * r * time;
  double eta = 1 / N * time;
  switch (type) {
    case SMC: {
      mat3 << -rho,  rho,   0,
                 0, -eta, eta,
                 0,    0,   0;
      return mat3.exp();
    }
    case SMC1: {
      mat3 << -rho,    rho,   0,
               eta, -2*eta, eta,
                 0,      0,   0;
      return mat3.exp();
    }
    case CSC: {
      four_dt mat4;
      mat4 << -rho,                  rho,        0,   0,
               eta, -(2 * eta + rho / 2),  rho / 2, eta,
                 0,              4 * eta, -5 * eta, eta,
                 0,                    0,        0,   0;
      return mat4.exp();
     }
     default:
      throw std::runtime_error("Unknown transition matrix requested.");
    }
}

double Transition::getTransitionFromStoInterval(double r, double timeS, double fromTime, double toTime,
                                                TransitionType type) {
  return getCumulativeTransitionProbability(r, timeS, toTime, type) // toCum
         - getCumulativeTransitionProbability(r, timeS, fromTime, type); // fromCum
}

std::vector<double> Transition::expectedIntervalTimesPiecewise() {
  std::vector<double> expectedTimes;
  for (unsigned i = 0; i < mDiscretization.size() - 1; i++) {
    expectedTimes.push_back(
        expectedTimeFromStoT(mDiscretization[i], mDiscretization[i + 1]));
  }
  return expectedTimes;
}

double Transition::expectedTimeFromStoT(double timeS, double timeT) {
  // TODO: what happens when findIntervalForTime returns -1?
  unsigned indexFrom = findIntervalForTime(timeS);
  unsigned indexTo = findIntervalForTime(timeT);
  double expected = 0.;
  double rate = 0.;
  for (unsigned i = indexFrom; i < indexTo + 1; i++) {
    double time0 = std::max(timeS, mTimeVectorPlusInfinity[i]);
    double time1 = std::min(timeT, mTimeVectorPlusInfinity[i + 1]);
    double N = mSize[i];
    double T = time1 - time0;
    if (time0 == time1) {
      continue;
    }
    double expectedThisPiece = (time1 == std::numeric_limits<double>::infinity())
      ? std::exp((timeS - time0) / N) * (N - timeS + time0)
      : std::exp(timeS / N) * ((N - timeS + time0) / std::exp(time0 / N) - (N - timeS + time1) / std::exp(time1 / N));
    rate -= T / N;
    expected += expectedThisPiece;
  }
  // prob having coalesced = 1 - prob not having coalesced
  return expected / (1 - std::exp(rate)) + timeS;
}

double Transition::getCumulativeTransitionProbability(double r, double timeS, double timeT, TransitionType type) {
  mat_dt Omega;
  if (timeT < timeS) {
    Omega = computeTransitionPiecewiseUpToTimeT(r, timeT, type);
    return Omega(
        0,
        (type == CSC) ? 3 : 2
    );
  } else if (timeT == timeS) {
    Omega = computeTransitionPiecewiseUpToTimeT(r, timeS, type);
    return Omega(0, 0) + Omega(
        0,
        (type == CSC) ? 3 : 2
    );
  } else {
    Omega = computeTransitionPiecewiseUpToTimeT(r, timeS, type);
    double cumCoalFromStoT = cumulativeCoalesceFromStoT(timeS, timeT);
    if (type == CSC) {
      return Omega(0, 0) + cumCoalFromStoT * (Omega(0, 1) + Omega(0, 2)) + Omega(0, 3);
    } else {
      return Omega(0, 0) + cumCoalFromStoT * Omega(0, 1) + Omega(0, 2);
    }
  }
}

double Transition::cumulativeCoalesceFromStoT(double timeS, double timeT) {
    return 1 - getSizeInPiecewiseAtTimeT(timeT) * coalesceFromStoT(timeS, timeT);
}

double Transition::getSizeInPiecewiseAtTimeT(double timeT) {
    return mSize[findIntervalForTime(timeT)];
}

double Transition::coalesceFromStoT(double timeS, double timeT) {
  if (timeT == std::numeric_limits<double>::infinity()) return 0.;

  unsigned indexFrom = findIntervalForTime(timeS);
  unsigned indexTo = findIntervalForTime(timeT);
  double rate = 0;
  for (auto i = indexFrom; i <= indexTo; i++)
    rate += (std::max(timeS, mTime[i]) - std::min(
          timeT, mTime[i + 1])) / mSize[i];
  double Nt = getSizeInPiecewiseAtTimeT(timeT);
  return 1 / Nt * std::exp(rate);
}

mat_dt Transition::computeTransitionPiecewiseUpToTimeT(double r, double time, TransitionType type) {
  unsigned indexTo = findIntervalForTime(time);
  mat_dt matrix = identity(type);
  for (unsigned i = 0; i <= indexTo - 1; i++) {
    matrix *= getExponentiatedTransitionMatrix(mSize[i], r, mTime[i + 1] - mTime[i], type);
  }
  matrix *= getExponentiatedTransitionMatrix(mSize[indexTo], r, time - mTime[indexTo], type);
  return matrix;
}

unsigned int Transition::findIntervalForTime(double t) {
  if (t == std::numeric_limits<double>::infinity()) return static_cast<unsigned int>(mSize.size()) - 1;

  for(unsigned i = 0; i < mSize.size(); i++) {
    if (t >= mTime[i] && t < mTime[i + 1]) return i;
  }
  throw std::runtime_error("Could not find interval for time: " + std::to_string(t));
}

double Transition::notCoalesceFromStoT(double timeS, double timeT) {
  if (timeT == std::numeric_limits<double>::infinity()) return 0.;
  unsigned indexFrom = findIntervalForTime(timeS);
  unsigned indexTo = findIntervalForTime(timeT);
  double rate = 0;
  for (unsigned i = indexFrom; i <= indexTo; i++) {
    rate += (std::max(timeS, mTime[i]) - std::min(timeT, mTime[i + 1])) / mSize[i];
  }
  return std::exp(rate);
}

mat_dt Transition::computeTransitionPiecewiseFromTimeSToTimeT(double r, double timeS, double timeT, TransitionType type) {
  mat_dt matrix = identity(type);
  unsigned indexFrom = findIntervalForTime(timeS);
  unsigned indexTo = findIntervalForTime(timeT);
  for (unsigned i = indexFrom; i <= indexTo; i++) {
    matrix *= getExponentiatedTransitionMatrix(mSize[i], r, (std::min(timeT, mTime[i + 1]) -
          std::max(timeS, mTime[i])), type);
  }
  return matrix;
}

double Transition::cumulativeCoalesceFromStoTsmart(double timeS, double timeT) {
  return 1 - notCoalesceFromStoT(timeS, timeT);
}

std::pair<mat_dt, mat_dt> Transition::getOmegas(double r, TransitionType type) {
  int cols = (type == CSC) ? 4 : 3;
  mat_dt omegasAtBoundaries, omegasAtExpectedTimes, latestOmega, M;
  latestOmega = identity(type);
  omegasAtBoundaries.resize(mStates + 1, cols);
  omegasAtExpectedTimes.resize(mStates, cols);
  omegasAtBoundaries.row(0) = latestOmega.row(0);
  for (unsigned i = 0; i < mDiscretization.size() - 1; i++) {
    double intervalStartTime = mDiscretization[i];
    double intervalExpTime = mExpectedTimes[i];
    double intervalEndTime = mDiscretization[i + 1];
    M = computeTransitionPiecewiseFromTimeSToTimeT(r, intervalStartTime, intervalExpTime, type);
    latestOmega *= M;
    omegasAtExpectedTimes.row(i) = latestOmega.row(0);
    if (intervalEndTime == std::numeric_limits<double>::infinity())
      M = identity(type);
    else
      M = computeTransitionPiecewiseFromTimeSToTimeT(r, intervalExpTime, intervalEndTime, type);
    latestOmega *= M;
    omegasAtBoundaries.row(i + 1) = latestOmega.row(0);
  }
  return std::make_pair(omegasAtBoundaries, omegasAtExpectedTimes);
}

void Transition::computeCoalescentVectors() {
  for (unsigned i = 0; i < mExpectedTimes.size(); i++) {
    double timeFrom = mDiscretization.at(i);
    double timeTo = mDiscretization.at(i + 1);
    double expTimeFrom = mExpectedTimes[i];
    if (i < mExpectedTimes.size() - 1) {
      mProbNotCoalesceBetweenExpectedTimes.push_back(notCoalesceFromStoT(expTimeFrom, mExpectedTimes[i + 1]));
    }
    mProbNotCoalesceBetweenTimeIntervals.push_back(notCoalesceFromStoT(timeFrom, timeTo));
    mProbCoalesceBetweenExpectedTimesAndUpperLimit.push_back(
      cumulativeCoalesceFromStoTsmart(expTimeFrom, timeTo));
  }
  // do U and RR up to states - 2
  mColumnRatios.resize(mStates - 1);
  // TODO: starting from 1, is this intended?
  for (unsigned i = 1; i < mStates - 1; i++) {
    double thisCR = mProbNotCoalesceBetweenTimeIntervals[i] *
      (1 - mProbNotCoalesceBetweenTimeIntervals[i + 1]) /
      (1 - mProbNotCoalesceBetweenTimeIntervals[i]);
    mColumnRatios(i) = std::isnan(thisCR) ? 1. : thisCR;
  }
}

std::vector<double> Transition::getCoalDist() {
  std::vector<double> coalDist;
  double lastCoal = 0.;
  for (double discretization : mDiscretization) {
    double coal = cumulativeCoalesceFromStoT(0., discretization);
    coalDist.push_back(coal - lastCoal);
    lastCoal = coal;
  }
  return coalDist;
}

} // namespace asmc
