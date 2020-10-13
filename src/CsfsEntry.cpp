// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include "CsfsEntry.hpp"

#include "EigenTypes.hpp"

#include <fmt/core.h>
#include <fmt/ranges.h>
#include <fmt/ostream.h>

#include <stdexcept>

namespace asmc {

CSFSEntry::CSFSEntry(std::vector<double> timeVector, std::vector<double> sizeVector, double mu, double from, double to, unsigned int samples,
                     mat_dt csfs)
    : mMu(mu), mFrom(from), mTo(to), mSamples(samples),
      mTimeVector(std::move(timeVector)),
      mSizeVector(std::move(sizeVector)), mCSFS(std::move(csfs)) {

  if (mTimeVector.size() != mSizeVector.size() || mCSFS.rows() != 3
      || mCSFS.cols() != mSamples - 1 || mFrom >= mTo) {

    const std::string errorMessage =
        fmt::format("Time vector:\n{}\nSize vector:\n{}\n{}\nFrom: {} To: {}\nMalformed CSFS entry.",
                    mTimeVector, mSizeVector, mCSFS, mFrom, mTo);
    throw std::runtime_error(errorMessage);
  }
}

std::string CSFSEntry::toString() const {
  return fmt::format("Time:\t{}\nSize:\t{}\nMu:\t{}\nSamples:\t{}\nInterval:\t{}\t{}\n{}\n",
      fmt::join(mTimeVector, " "), fmt::join(mSizeVector, " "), mMu, mSamples, mFrom, mTo, mCSFS);
}

} // namespace asmc
