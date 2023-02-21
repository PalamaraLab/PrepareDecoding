// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include "DefaultDemographies.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <vector>

namespace fs = std::filesystem;

using Catch::Matchers::WithinAbs;

namespace asmc {

TEST_CASE("Default demographies can be written to file", "[DefaultDemographies]") {

  // Check CEU by writing it to file and checking against the built-in
  {
    demo::saveDemography(PREPARE_DECODING_TEST_DIR "/data", "CEU");

    auto [times, sizes] = demo::getBuiltInDemography("CEU");

    std::ifstream file(PREPARE_DECODING_TEST_DIR "/data/CEU.demo");
    std::string line;
    double time{};
    double size{};
    for (auto i = 0ul; i < times.size(); ++i) {
      std::getline(file, line);
      std::stringstream ss(line);
      ss >> time >> size;
      CHECK_THAT(times.at(i), WithinAbs(time, 1e-12));
      CHECK_THAT(sizes.at(i), WithinAbs(size, 1e-12));
    }

    std::getline(file, line);
    CHECK(line.empty());

    fs::remove(PREPARE_DECODING_TEST_DIR "/data/CEU.demo");
  }

}

} // namespace asmc
