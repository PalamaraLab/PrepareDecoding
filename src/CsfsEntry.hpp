// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#ifndef PREPAREDECODING_CSFSENTRY_HPP
#define PREPAREDECODING_CSFSENTRY_HPP

#include "EigenTypes.hpp"

namespace asmc {

class CsfsEntry {

private:
  double mMu = {};
  double mFrom = {};
  double mTo = {};

  int mSamples = {};

  array_dt mTimeVector = {};
  array_dt mSizeVector = {};

  mat_dt mCsfs = {};

public:
  CsfsEntry(array_dt timeVector, array_dt sizeVector, double mu, double from, double to, int samples, mat_dt csfs);

  std::string toString();
};

} // namespace asmc

#endif // PREPAREDECODING_CSFSENTRY_HPP
