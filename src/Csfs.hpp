// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#ifndef PREPAREDECODING_CSFS_HPP
#define PREPAREDECODING_CSFS_HPP

#include <map>
#include "EigenTypes.hpp"
#include "Transition.hpp"
#include "CsfsEntry.hpp"
#include "ArraySpectrum.hpp"

namespace asmc {

enum class CSFSParserState { Null, Time, Size, Mu, Samples, Interval, CSFS };
const std::map<std::string, CSFSParserState> stateMap = {
  {"Size:", CSFSParserState::Size},
  {"Time:", CSFSParserState::Time},
  {"Mu:", CSFSParserState::Mu},
  {"Samples:", CSFSParserState::Samples},
  {"Interval:", CSFSParserState::Interval}
};

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
  unsigned int mSamples = 0;

public:

  CSFS() = default;
  explicit CSFS(std::map<double, CSFSEntry> CSFS_);
  static CSFSParserState currentState(const std::string& line);
  static std::pair<CSFSParserState, int> nextState(CSFSParserState state, int line);
  static CSFS loadFromFile(std::string_view filename);
  bool verify(std::vector<double> timeVectorOriginal, std::vector<double> sizeVectorOriginal,
      double mu, unsigned int samples, std::vector<double> discretizationOriginal);
  std::string toString() const;
  void fixAscertainment(Data data, unsigned int samples, Transition transition);
  static mat_dt computeClassicEmission(std::vector<double> expectedTimes, double mu);
  void computeArraySamplingFactors(Data data, unsigned int samples, Transition transition);
  void applyFactors();
  static std::map<double, CSFSEntry> foldCSFS(std::map<double, CSFSEntry> csfsMap);
  static std::pair<unsigned int, unsigned int> getFoldedObservationFromUnfolded(
      std::pair<unsigned int, unsigned int> unfolded, unsigned int totalSamples);
  static mat_dt compressCSFS(std::map<double, CSFSEntry> csfsMap);

  std::map<double, CSFSEntry>& getCSFS() { return mCSFS; }
  std::map<double, CSFSEntry>& getFoldedCSFS() { return mFoldedCSFS; }
  std::map<double, CSFSEntry>& getAscertainedCSFS() { return mAscertainedCSFS; }
  std::map<double, CSFSEntry>& getFoldedAscertainedCSFS() { return mFoldedAscertainedCSFS; }
  mat_dt& getCompressedAscertainedEmissionTable() { return mCompressedAscertainedEmissionTable; }
  unsigned int getSamples() const { return mSamples; }

};

} // namespace asmc

#endif // PREPAREDECODING_CSFS_HPP
