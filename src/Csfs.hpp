// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#ifndef PREPAREDECODING_CSFS_HPP
#define PREPAREDECODING_CSFS_HPP

#include "EigenTypes.hpp"
#include "ArraySpectrum.hpp"

namespace asmc {

enum class CSFSParserState { Null, Time, Size, Mu, Samples, Interval, CSFS };

class CSFS {

private:

  ArraySpectrum mArraySpectrum;
  std::map<double, CSFSEntry> mCSFS;
  std::map<double, CSFSEntry> mFoldedCSFS;

  // ascertained emission
  std::map<double, CSFSEntry> mAscertainedCSFS;
  std::map<double, CSFSEntry> mFoldedAscertainedCSFS;
  mat_dt mCompressedAscertainedEmissionTable;

  array_dt mArraySamplingFactors;
  int mSamples;
  int mCSFSLine = 0;

public:

  explicit CSFS(std::map<double, CSFSEntry> CSFS_);
  const static std::map<std::string, CSFSParserState> stateMap = {
    {"Size:", CSFSParserState::Size},
    {"Time:", CSFSParserState::Time},
    {"Mu:", CSFSParserState::Mu},
    {"Samples:", CSFSParserState::Samples},
    {"Interval:", CSFSParserState::Interval}
  };
  static CSFSParserState currentState(std::string line);
  static CSFSParserState nextState(CSFSParserState state);
  static CSFS loadFromFile(std::string_view filename);
  bool verify(std::vector<double> timeVectorOriginal, std::vector<double> sizeVectorOriginal,
      double mu, int samples, std::vector<double> discretizationOriginal);
  std::string toString();
  void fixAscertainment(Data data, int samples, Transition transition);
  static mat_dt computeClassicEmission(std::vector<double> expectedTimes, double mu);
  void computeArraySamplingFactors(Data data, int samples, Transition transition);
  void applyFactors();
  std::map<double, CSFSEntry> foldCSFS(std::map<double, CSFSEntry> CSFS);
  static std::pair<int, int> getFoldedObservationFromUnfolded(
      std::pair<int, int> unfolded, int totalSamples);
  mat_dt compressCSFS(std::map<double, CSFSEntry> CSFS);

};

} // namespace asmc

#endif // PREPAREDECODING_CSFS_HPP
