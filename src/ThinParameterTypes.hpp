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
} // namespace asmc

#endif // PREPAREDECODING_THIN_PARAMETER_TYPES_HPP
