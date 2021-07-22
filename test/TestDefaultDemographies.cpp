// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include "DefaultDemographies.hpp"

#include <catch2/catch.hpp>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <vector>

namespace fs = std::filesystem;

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
      CHECK(times.at(i) == Approx(time).epsilon(1e-12));
      CHECK(sizes.at(i) == Approx(size).epsilon(1e-12));
    }

    std::getline(file, line);
    CHECK(line.empty());

    fs::remove(PREPARE_DECODING_TEST_DIR "/data/CEU.demo");
  }

}

} // namespace asmc
