// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#ifndef PREPAREDECODING_CSFSENTRY_HPP
#define PREPAREDECODING_CSFSENTRY_HPP

#include <Eigen/Core>
#include <fmt/core.h>
#include <fmt/ostream.h>

#include <stdexcept>

namespace asmc {

using array_t = Eigen::ArrayXd;
using mat_t = Eigen::MatrixXd;

class CsfsEntry {

private:
  double mMu = {};
  double mFrom = {};
  double mTo = {};

  int mSamples = {};

  array_t mTimeVector = {};
  array_t mSizeVector = {};

  mat_t mCsfs = {};

public:
  CsfsEntry(array_t timeVector, array_t sizeVector, double mu, double from, double to, int samples, mat_t csfs)
      : mMu(mu), mFrom(from), mTo(to), mSamples(samples), mTimeVector(std::move(timeVector)),
        mSizeVector(std::move(sizeVector)), mCsfs(std::move(csfs)) {

    if (mTimeVector.size() != mSizeVector.size() || mCsfs.rows() != 3 || mCsfs.cols() != mSamples - 1 || mFrom >= mTo) {

      const std::string errorMessage =
          fmt::format("Time vector:\n{}\nSize vector:\n{}\nCSFS:\n{}\nFrom: {} To: {}\nMalformed CSFS entry.",
                      mTimeVector.transpose(), mSizeVector.transpose(), mCsfs, mFrom, mTo);
      throw std::runtime_error(errorMessage);
    }
  }

  std::string toString() {
    return fmt::format("Time:\t{}\nSize:\t{}\nMu:\t{}\nSamples:\t{}\nInterval:\t{}\t{}\n{}\n", mTimeVector.transpose(),
                       mSizeVector.transpose(), mMu, mSamples, mFrom, mTo, mCsfs);
  }
};

} // namespace asmc

#endif // PREPAREDECODING_CSFSENTRY_HPP
