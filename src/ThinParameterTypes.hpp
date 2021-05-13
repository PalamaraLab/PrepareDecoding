// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#ifndef PREPAREDECODING_THIN_PARAMETER_TYPES_HPP
#define PREPAREDECODING_THIN_PARAMETER_TYPES_HPP

#include <fmt/core.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <array>
#include <filesystem>
#include <string>
#include <string_view>

namespace fs = std::filesystem;

namespace asmc {

class Demography {
private:
  bool mFile = false;
  bool mBuiltIn = false;
  std::string mDemography;

  constexpr static std::array mValidDefaults{"ACB", "ASW", "BEB", "CDX", "CEU", "CHB", "CHS", "CLM", "ESN",
                                             "FIN", "GBR", "GIH", "GWD", "IBS", "ITU", "JPT", "KHV", "LWK",
                                             "MSL", "MXL", "PEL", "PJL", "PUR", "STU", "TSI", "YRI"};

public:
  Demography() : mBuiltIn{true}, mDemography{"CEU"} {
    fmt::print("Default demography set to CEU\n");
  }

  explicit Demography(std::string_view demography) : mDemography{demography} {

    if (fs::exists(demography) && fs::is_regular_file(demography)) {
      mFile = true;
      fmt::print("Demography set to file {}\n", mDemography);
    } else if (std::find(std::begin(mValidDefaults), std::end(mValidDefaults), demography) !=
               std::end(mValidDefaults)) {
      mBuiltIn = true;
      fmt::print("Demography set to built-in {}\n", mDemography);
    } else {
      auto error_message = fmt::format("Expected either a valid demography file, or a default (one of {}), but got {}",
                                       mValidDefaults, mDemography);
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

} // namespace asmc

#endif // PREPAREDECODING_THIN_PARAMETER_TYPES_HPP
