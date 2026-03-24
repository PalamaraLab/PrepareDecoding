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

#include "Transition.hpp"
#include "EigenTypes.hpp"

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include <unsupported/Eigen/MatrixFunctions>

namespace asmc {

std::string getTransitionTypeString(TransitionType tt) {
  switch (tt) {
  case TransitionType::SMC:
    return "SMC";
  case TransitionType::SMC1:
    return "SMC1";
  case TransitionType::CSC:
    return "CSC";
  }
  return "";
}

std::vector<double> Transition::getTimeExponentialQuantiles(int numQuantiles, const std::vector<double>& timeVector,
                                                            const std::vector<double>& sizeFromVector) {
  assert(!timeVector.empty());
  assert(timeVector.size() == sizeFromVector.size());
  assert(timeVector.back() != std::numeric_limits<double>::infinity());

  std::vector<double> timesWithInf = timeVector;
  timesWithInf.push_back(std::numeric_limits<double>::infinity());

  const double timeStep = 0.1;

  double pNotCoal = 1.0;
  double nextQuantile = 1.0 / numQuantiles;
  std::vector<double> quantiles{0.0};

  for (auto i = 0ul; i < timesWithInf.size() - 1; i++) {
    const double tStart = timesWithInf.at(i);
    const double tEnd = timesWithInf.at(i + 1);
    const double notCoalRate = 1.0 - timeStep / sizeFromVector.at(i);

    unsigned int count = 0u;
    double t = tStart;
    while (t < tEnd) {
      double newPNotCoal = pNotCoal * std::pow(notCoalRate, count);

      if (1.0 - newPNotCoal > nextQuantile) {
        quantiles.push_back(t);
        nextQuantile = static_cast<double>(quantiles.size()) / static_cast<double>(numQuantiles);

        if (quantiles.size() == static_cast<std::size_t>(numQuantiles)) {
          return quantiles;
        }
      }
      count++;
      t = tStart + count * timeStep;
    }

    pNotCoal *= std::pow(notCoalRate, count - 1);
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
  const double MAX_T = sizeFromVector.back() * 20;
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
  mTimeVectorPlusInfinity = mTime;
  mTimeVectorPlusInfinity.push_back(std::numeric_limits<double>::infinity());
  mExpectedTimes = expectedIntervalTimesPiecewise();
  mStates = static_cast<unsigned int>(mDiscretization.size()) - 1;
  computeCoalescentVectors();
}

mat_dt Transition::identity(TransitionType type) {
  if (type == TransitionType::CSC) {
    return four_dt::Identity();
  } else {
    return three_dt::Identity();
  }
}

unsigned int Transition::findIntervalForTime(double t) {
  if (t == std::numeric_limits<double>::infinity())
    return static_cast<unsigned int>(mSize.size()) - 1;
  auto it = std::upper_bound(mTime.begin(), mTime.end(), t);
  if (it == mTime.begin())
    throw std::runtime_error("Could not find interval for time: " + std::to_string(t));
  return static_cast<unsigned int>(std::distance(mTime.begin(), it)) - 1;
}

mat_dt Transition::getExponentiatedTransitionMatrix(double N, double r, double time, TransitionType type) {
  double rho = 2.0 * r * time;
  double eta = time / N;
  switch (type) {
    case TransitionType::SMC: {
      // Upper triangular: closed-form exponential
      double er = std::exp(-rho);
      double en = std::exp(-eta);
      double m01;
      if (rho == 0.0) {
        m01 = 0.0;
      } else if (std::abs(eta - rho) > 1e-10 * std::max(rho, eta)) {
        m01 = rho * (er - en) / (eta - rho);
      } else {
        m01 = rho * er;
      }
      three_dt result;
      result << er,  m01,         1.0 - er - m01,
                0.0, en,          1.0 - en,
                0.0, 0.0,         1.0;
      return result;
    }
    case TransitionType::SMC1: {
      // 2x2 eigendecomposition for B = [[-rho, rho], [eta, -2*eta]], row-sum for 3rd column
      if (rho == 0.0 && eta == 0.0) return three_dt::Identity();
      double s = std::sqrt(rho * rho + 4.0 * eta * eta);
      double lam1 = (-(rho + 2.0 * eta) + s) / 2.0;
      double lam2 = (-(rho + 2.0 * eta) - s) / 2.0;
      double e1 = std::exp(lam1);
      double e2 = std::exp(lam2);
      // exp(B) = (e1*(B - lam2*I) - e2*(B - lam1*I)) / (lam1 - lam2)
      double b00 = (e1 * (-rho - lam2) - e2 * (-rho - lam1)) / s;
      double b01 = rho * (e1 - e2) / s;
      double b10 = eta * (e1 - e2) / s;
      double b11 = (e1 * (-2.0 * eta - lam2) - e2 * (-2.0 * eta - lam1)) / s;
      three_dt result;
      result << b00, b01, 1.0 - b00 - b01,
                b10, b11, 1.0 - b10 - b11,
                0.0, 0.0, 1.0;
      return result;
    }
    case TransitionType::CSC: {
      // 3x3 block + absorbing state: compute exp(B) for non-absorbing block,
      // derive 4th column from row-sum = 1 (exact, avoids 4x4 exponential)
      three_dt B;
      B << -rho,                     rho,           0.0,
            eta,  -(2.0 * eta + rho / 2.0),  rho / 2.0,
            0.0,               4.0 * eta,    -5.0 * eta;
      three_dt expB = B.exp();
      four_dt result;
      result << expB(0,0), expB(0,1), expB(0,2), 1.0 - expB.row(0).sum(),
                expB(1,0), expB(1,1), expB(1,2), 1.0 - expB.row(1).sum(),
                expB(2,0), expB(2,1), expB(2,2), 1.0 - expB.row(2).sum(),
                0.0,       0.0,       0.0,       1.0;
      return result;
    }
    default:
      throw std::runtime_error("Unknown transition matrix requested.");
  }
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
  return expected / (1 - std::exp(rate)) + timeS;
}

double Transition::getSizeInPiecewiseAtTimeT(double timeT) {
  return mSize[findIntervalForTime(timeT)];
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

double Transition::cumulativeCoalesceFromStoT(double timeS, double timeT) {
  if (timeT == std::numeric_limits<double>::infinity()) return 1.0;
  unsigned indexFrom = findIntervalForTime(timeS);
  unsigned indexTo = findIntervalForTime(timeT);
  double rate = 0;
  for (auto i = indexFrom; i <= indexTo; i++)
    rate += (std::max(timeS, mTime[i]) - std::min(timeT, mTime[i + 1])) / mSize[i];
  double Nt = mSize[indexTo];
  return 1.0 - Nt * (1.0 / Nt * std::exp(rate));
}

mat_dt Transition::computeTransitionPiecewiseFromTimeSToTimeT(double r, double timeS, double timeT, TransitionType type) {
  mat_dt matrix = identity(type);
  unsigned indexFrom = findIntervalForTime(timeS);
  unsigned indexTo = findIntervalForTime(timeT);
  for (unsigned i = indexFrom; i <= indexTo; i++) {
    matrix *= getExponentiatedTransitionMatrix(mSize[i], r,
        std::min(timeT, mTime[i + 1]) - std::max(timeS, mTime[i]), type);
  }
  return matrix;
}

std::pair<mat_dt, mat_dt> Transition::getOmegas(double r, TransitionType type) {
  int cols = (type == TransitionType::CSC) ? 4 : 3;
  mat_dt omegasAtBoundaries(mStates + 1, cols);
  mat_dt omegasAtExpectedTimes(mStates, cols);

  mat_dt latestOmega = identity(type);
  omegasAtBoundaries.row(0) = latestOmega.row(0);

  for (unsigned i = 0; i < mStates; i++) {
    latestOmega *= computeTransitionPiecewiseFromTimeSToTimeT(r, mDiscretization[i], mExpectedTimes[i], type);
    omegasAtExpectedTimes.row(i) = latestOmega.row(0);
    if (mDiscretization[i + 1] != std::numeric_limits<double>::infinity()) {
      latestOmega *= computeTransitionPiecewiseFromTimeSToTimeT(r, mExpectedTimes[i], mDiscretization[i + 1], type);
    }
    omegasAtBoundaries.row(i + 1) = latestOmega.row(0);
  }

  return {omegasAtBoundaries, omegasAtExpectedTimes};
}

std::tuple<vec_dt, vec_dt, vec_dt, vec_dt>
Transition::getLinearTimeDecodingQuantitiesAndMatrixGivenDistance(double rho) {
  auto omegas = getOmegas(rho, mType);
  mat_dt& omegasAtBoundaries = omegas.first;
  mat_dt& omegasAtExpectedTimes = omegas.second;
  vec_dt D(mStates), B(mStates - 1), U(mStates - 1), RR(mStates - 1);
  D.setZero(); B.setZero(); U.setZero(); RR.setZero();

  const int lastCol = (mType == TransitionType::CSC) ? 3 : 2;
  const bool isCSC = (mType == TransitionType::CSC);

  auto omegaS = [&](unsigned i) -> double {
    return isCSC ? omegasAtExpectedTimes(i, 1) + omegasAtExpectedTimes(i, 2)
                 : omegasAtExpectedTimes(i, 1);
  };

  for (unsigned i = 0; i < mStates; i++) {
    D[i] = omegasAtExpectedTimes(i, 0)
         + mProbCoalesceBetweenExpectedTimesAndUpperLimit[i] * omegaS(i)
         + omegasAtExpectedTimes(i, lastCol) - omegasAtBoundaries(i, lastCol);
    if (i > 0)
      B[i - 1] = omegasAtBoundaries(i, lastCol) - omegasAtBoundaries(i - 1, lastCol);
  }

  for (unsigned i = 0; i < mStates - 1; i++) {
    double oS = omegaS(i);
    U[i] = oS * (1.0 - mProbCoalesceBetweenExpectedTimesAndUpperLimit[i])
             * (1.0 - mProbNotCoalesceBetweenTimeIntervals[i + 1]);
    if (i < mStates - 2)
      RR[i] = (rho == 0.0) ? 1.0 : oS * mProbNotCoalesceBetweenExpectedTimes[i] / omegaS(i + 1);
  }

  return {D, B, U, RR};
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
        1.0 - notCoalesceFromStoT(expTimeFrom, timeTo));
  }
  mColumnRatios.resize(mStates - 1);
  for (unsigned i = 1; i < mStates - 1; i++) {
    double thisCR = mProbNotCoalesceBetweenTimeIntervals[i] *
      (1.0 - mProbNotCoalesceBetweenTimeIntervals[i + 1]) /
      (1.0 - mProbNotCoalesceBetweenTimeIntervals[i]);
    mColumnRatios(i) = std::isnan(thisCR) ? 1.0 : thisCR;
  }
}

std::vector<double> Transition::getCoalDist() {
  std::vector<double> coalDist;
  double lastCoal = 0.;
  for (unsigned i = 1; i < mDiscretization.size(); i++) {
    double coal = cumulativeCoalesceFromStoT(0., mDiscretization[i]);
    coalDist.push_back(coal - lastCoal);
    lastCoal = coal;
  }
  return coalDist;
}

} // namespace asmc
