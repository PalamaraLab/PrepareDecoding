// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include "CsfsEntry.hpp"

#include "EigenTypes.hpp"

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <stdexcept>

namespace asmc {

CSFSEntry::CSFSEntry(array_dt timeVector, array_dt sizeVector, double mu, double from, double to, int samples,
                     mat_dt csfs)
    : mMu(mu), mFrom(from), mTo(to), mSamples(samples),
      mTimeVector(std::move(timeVector)),
      mSizeVector(std::move(sizeVector)), mCSFS(std::move(csfs)) {

  if (mTimeVector.size() != mSizeVector.size() || mCSFS.rows() != 3
      || mCSFS.cols() != mSamples - 1 || mFrom >= mTo) {

    const std::string errorMessage =
        fmt::format("Time vector:\n{}\nSize vector:\n{}\nCSFS:\n{}\nFrom: {} To: {}\nMalformed CSFS entry.",
                    mTimeVector.transpose(), mSizeVector.transpose(), mCSFS, mFrom, mTo);
    throw std::runtime_error(errorMessage);
  }
}

std::string CSFSEntry::toString() {
  return fmt::format("Time:\t{}\nSize:\t{}\nMu:\t{}\nSamples:\t{}\nInterval:\t{}\t{}\n{}\n", mTimeVector.transpose(),
                     mSizeVector.transpose(), mMu, mSamples, mFrom, mTo, mCSFS);
}

} // namespace asmc
