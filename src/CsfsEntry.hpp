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

  unsigned int mSamples = {};

  std::vector<double> mTimeVector = {};
  std::vector<double> mSizeVector = {};

  mat_dt mCSFS = {};

public:
  CSFSEntry(std::vector<double> timeVector, std::vector<double> sizeVector,
      double mu, double from, double to, unsigned int samples, mat_dt csfs);

  std::vector<double>& getTime() { return mTimeVector; }
  std::vector<double>& getSize() { return mSizeVector; }
  double getMu() { return mMu; }
  double getFrom() { return mFrom; }
  double getTo() { return mTo; }
  unsigned int getSamples() { return mSamples; }
  const mat_dt& getCSFSMatrix() const { return mCSFS; }
  void setCSFSMatrix(const mat_dt& csfs) { mCSFS = csfs; }

  std::string toString() const;
};

} // namespace asmc

#endif // PREPAREDECODING_CSFSENTRY_HPP
