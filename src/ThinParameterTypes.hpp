// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#ifndef PREPAREDECODING_THIN_PARAMETER_TYPES_HPP
#define PREPAREDECODING_THIN_PARAMETER_TYPES_HPP

#include "DefaultDemographies.hpp"

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <array>
#include <filesystem>
#include <string>
#include <string_view>
#include <utility>

namespace fs = std::filesystem;

namespace asmc {

class Demography {
private:
  bool mFile = false;
  bool mBuiltIn = false;
  std::string mDemography;

public:
  Demography() : mBuiltIn{true}, mDemography{"CEU"} {
    fmt::print("Default demography set to CEU\n");
  }

  explicit Demography(std::string_view demography) : mDemography{demography} {

    if (fs::exists(demography) && fs::is_regular_file(demography)) {
      mFile = true;
      fmt::print("Demography set to file {}\n", mDemography);
    } else if (demo::isValidDemography(mDemography)) {
      mBuiltIn = true;
      fmt::print("Demography set to built-in {}\n", mDemography);
    } else {
      auto error_message = fmt::format("Expected either a valid demography file, or a default (one of {}), but got {}",
                                       demo::validDemographies, mDemography);
      throw std::runtime_error(error_message);
    }
  }

  [[nodiscard]] bool isFile() const {
    return mFile;
  }
  [[nodiscard]] bool isBuiltIn() const {
    return mBuiltIn;
  }
  [[nodiscard]] const std::string& getDemography() const {
    return mDemography;
  }
};

class Discretization {
private:
  bool mFile = false;
  std::string mDiscretizationFile;
  std::vector<double> mDiscretizationPoints;
  int mNumAdditionalPoints = 0;

public:
  explicit Discretization(std::string_view discretizationFile) : mDiscretizationFile{discretizationFile} {

    if (fs::exists(discretizationFile) && fs::is_regular_file(discretizationFile)) {
      mFile = true;
      fmt::print("Discretization file set to file {}\n", mDiscretizationFile);
    } else {
      auto error_message = fmt::format("Expected a valid discretization file, but got {}", mDiscretizationFile);
      throw std::runtime_error(error_message);
    }
  }

  Discretization(std::vector<double> discretizationPoints, const int numAdditionalPoints)
      : mDiscretizationPoints{std::move(discretizationPoints)}, mNumAdditionalPoints{numAdditionalPoints} {
    bool valid = true;
    if (mDiscretizationPoints.empty() && mNumAdditionalPoints < 2) {
      valid = false;
    } else if (!mDiscretizationPoints.empty() && mDiscretizationPoints.front() != 0.0) {
      valid = false;

      for (auto i = 1ul; i < mDiscretizationPoints.size(); ++i) {
        if (mDiscretizationPoints.at(i) <= mDiscretizationPoints.at(i - 1)) {
          valid = false;
        }
      }
    }

    if (!valid) {
      auto error_message = fmt::format(
          "Expected a monotonic increasing vector of discretization points, starting with 0.0, and a number "
          "of additional points to calculate, but got {} and {}",
          mDiscretizationPoints, mNumAdditionalPoints);
      throw std::runtime_error(error_message);
    }
  }

  Discretization(const std::vector<std::pair<double, int>>& discretizationPointPairs, const int numAdditionalPoints)
      : mNumAdditionalPoints{numAdditionalPoints} {

    // Check for sensible inputs
    bool valid = true;
    if (discretizationPointPairs.empty() && numAdditionalPoints < 2) {
      valid = false;
    }

    if (!discretizationPointPairs.empty()) {
      mDiscretizationPoints = {0.0};
      for (const auto& [val, num] : discretizationPointPairs) {

        // Check for monotonicity
        if (val < 0.0 || num < 1) {
          valid = false;
          break;
        }

        const double baseValue = mDiscretizationPoints.back();
        for (int i = 0; i < num; ++i) {
          mDiscretizationPoints.push_back(baseValue + (1 + i) * val);
        }
      }
    }

    if (!valid) {
      auto error_message =
          fmt::format("Expected pairs of form [val (double), num (int)] of discretization points and a number "
                      "of additional points to calculate, but got {} and {}",
                      discretizationPointPairs, numAdditionalPoints);
      throw std::runtime_error(error_message);
    }
  }

  [[nodiscard]] bool isFile() const {
    return mFile;
  }

  [[nodiscard]] const std::string& getDiscretizationFile() const {
    return mDiscretizationFile;
  }

  [[nodiscard]] const std::vector<double>& getDiscretizationPoints() const {
    return mDiscretizationPoints;
  }

  [[nodiscard]] int getNumAdditionalPoints() const {
    return mNumAdditionalPoints;
  }
};

class Frequencies {
private:
  bool mFile = false;
  bool mBuiltIn = false;
  bool mGz = false;
  std::string mFreqIdentifier;
  unsigned mNumSamples = 0;

  constexpr static std::array mValidDefaults{"UKBB"};
  constexpr static std::array mValidSamples{50, 100, 200, 300};

public:

  Frequencies() = default;

  explicit Frequencies(std::string_view fileRoot) : mFreqIdentifier{fileRoot} {

    // Check whether we're passed a file root containing a frequencies file
    if (std::string freqGzFile = fmt::format("{}.frq.gz", mFreqIdentifier);
        fs::exists(freqGzFile) && fs::is_regular_file(freqGzFile)) {
      mFile = true;
      mGz = true;
      mFreqIdentifier = freqGzFile;
      fmt::print("Frequencies file set to {}\n", mFreqIdentifier);
      return;
    }

    if (std::string freqFile = fmt::format("{}.frq", mFreqIdentifier);
        fs::exists(freqFile) && fs::is_regular_file(freqFile)) {
      mFile = true;
      mFreqIdentifier = freqFile;
      fmt::print("Frequencies file set to {}\n", mFreqIdentifier);
      return;
    }

  }

  Frequencies(std::string_view frequencies, const unsigned numSamples)
      : mFreqIdentifier{frequencies}, mNumSamples{numSamples} {

    // First check if it's a valid built-in
    bool validDataSource =
        std::find(std::begin(mValidDefaults), std::end(mValidDefaults), mFreqIdentifier) != std::end(mValidDefaults);
    bool validNumSamples =
        std::find(std::begin(mValidSamples), std::end(mValidSamples), numSamples) != std::end(mValidSamples);

    if (validDataSource && validNumSamples) {
      mBuiltIn = true;
      fmt::print("Frequencies set to built-in {} with {} samples\n", mFreqIdentifier, numSamples);
      return;
    }

    // If we're passed a valid file, use that
    if (fs::exists(mFreqIdentifier) && fs::is_regular_file(mFreqIdentifier)) {
      mFile = true;
      fmt::print("Frequencies file set to {}\n", mFreqIdentifier);
      return;
    }

    // If neither of the above, it's an error
    auto error_message = fmt::format("Expected either a valid frequencies file, or a valid data source (one of {}) and "
                                     "a valid number of samples (one of {}), but got {} with {} samples",
                                     mValidDefaults, mValidSamples, mFreqIdentifier, mNumSamples);
    throw std::runtime_error(error_message);
  }

  [[nodiscard]] bool isFile() const {
    return mFile;
  }

  [[nodiscard]] bool isBuiltIn() const {
    return mBuiltIn;
  }

  [[nodiscard]] bool isGz() const {
    return mGz;
  }

  [[nodiscard]] const std::string& getFreqIdentifier() const {
    return mFreqIdentifier;
  }

  [[nodiscard]] unsigned getNumSamples() const {
    return mNumSamples;
  }
};

} // namespace asmc

#endif // PREPAREDECODING_THIN_PARAMETER_TYPES_HPP
