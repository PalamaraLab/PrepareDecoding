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
#include <unsupported/Eigen/MatrixFunctions>
#include "EigenTypes.hpp"
#include "Transition.hpp"

namespace asmc {
    
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
  mTimeVectorPlusInfinity = timeVector;
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
  return getCumulativeTransitionPobability(r, timeS, toTime, type) // toCum
         - getCumulativeTransitionPobability(r, timeS, fromTime, type); // fromCum
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

double Transition::getCumulativeTransitionPobability(double r, double timeS, double timeT, TransitionType type) {
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
    double timeFrom = mDiscretization[i];
    double timeTo = mDiscretization[i + 1];
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
