// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#ifndef PREPAREDECODING_DATA_HPP
#define PREPAREDECODING_DATA_HPP

#include <vector>
#include "EigenTypes.hpp"

namespace asmc {

class Data {

private:

  std::vector<double> mAllSNPsFreq = {};
  std::vector<unsigned int> mAllSNPsMinorAlleles = {};
  std::vector<unsigned int> mAllSNPsAlleleCounts = {};
  unsigned int mHaploidSampleSize = 0;

  void readMinorAlleleFrequencies(std::string_view freqFile);
  void readMinorAlleleFrequenciesGz(std::string_view freqFile);
  void readMinorAlleleFrequenciesLine(const std::string& line);

  void computeMinorAlleleFrequenciesFromHaps(std::string_view hapsFileRoot);

  static std::string identifyAppropriateHapsFile(std::string_view hapsFileRoot);

public:

  /**
   * Default constructor taking no arguments.
   *
   * It is assumed that the use of this constructor is followed by use of the addFreq() method.
   */
  Data() = default;

  /**
   * Construct a Data object from the given file root.
   *
   * @param hapsFileRoot the root location of the .frq.gz, .frq, or haps file
   */
  explicit Data(std::string_view hapsFileRoot);

  unsigned int getHaploidSampleSize() const { return mHaploidSampleSize; }
  std::vector<double> getAllSNPsFreq();
  std::vector<unsigned int> getAllSNPsMinorAlleles();
  std::vector<unsigned int> getAllSNPsAlleleCounts();

  void addFreq(std::string_view freqFile);

};

} // namespace asmc

#endif // PREPAREDECODING_DATA_HPP
