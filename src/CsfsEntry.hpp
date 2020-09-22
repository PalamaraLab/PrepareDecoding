// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#ifndef PREPAREDECODING_CSFSENTRY_HPP
#define PREPAREDECODING_CSFSENTRY_HPP

#include "EigenTypes.hpp"

namespace asmc {

class CSFSEntry {

private:
  double mMu = {};
  double mFrom = {};
  double mTo = {};

  int mSamples = {};

  std::vector<double> mTimeVector = {};
  std::vector<double> mSizeVector = {};

  mat_dt mCSFS = {};

public:
  CSFSEntry(std::vector<double> timeVector, std::vector<double> sizeVector,
      double mu, double from, double to, int samples, mat_dt csfs);

  std::vector<double>& getTime() { return mTimeVector; }
  std::vector<double>& getSize() { return mSizeVector; }
  double getMu() { return mMu; }
  double getFrom() { return mFrom; }
  double getTo() { return mTo; }
  int getSamples() { return mSamples; }

  std::string toString();
};

} // namespace asmc

#endif // PREPAREDECODING_CSFSENTRY_HPP
