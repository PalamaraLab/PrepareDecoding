// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include "ThinParameterTypes.hpp"

#include <catch2/catch.hpp>

namespace asmc {

TEST_CASE("Test demography type", "[ThinParameterTypes]") {

  // Default
  {
    const Demography d;
    CHECK(d.isBuiltIn());
    CHECK(!d.isFile());
    CHECK(d.getDemography() == "CEU");
  }

  // Valid built in
  {
    const Demography d("ACB");
    CHECK(d.isBuiltIn());
    CHECK(!d.isFile());
    CHECK(d.getDemography() == "ACB");
  }

  // Invalid built in
  {
    CHECK_THROWS_WITH(Demography("ZZZ"), Catch::Contains("Expected either a valid demography"));
    CHECK_THROWS_WITH(Demography("ZZZ"), Catch::Contains("ZZZ"));
  }

  // Valid file
  {
    const Demography d(PREPARE_DECODING_TEST_DIR "/data/test.demo");
    CHECK(!d.isBuiltIn());
    CHECK(d.isFile());
    CHECK(Catch::contains(d.getDemography(), "test.demo"));
  }

  // Invalid file
  {
    CHECK_THROWS_WITH(Demography(PREPARE_DECODING_TEST_DIR "/data/test.doesNotExist"),
                      Catch::Contains("test.doesNotExist"));
  }
}

TEST_CASE("Test discretization type", "[ThinParameterTypes]") {

  // Construct with vector of values
  {
    Discretization d({0.0, 1.0, 2.0, 3.0}, 12);
    CHECK(d.getNumAdditionalPoints() == 12);
    CHECK(d.getDiscretizationPoints() == std::vector<double>{0.0, 1.0, 2.0, 3.0});
  }

  // Construct with vector of values not starting with zero
  {
    CHECK_THROWS_WITH(Discretization({1.0, 2.0, 3.0}, 12),
                      Catch::Contains("Expected a monotonic increasing vector of discretization points"));
  }

  // Construct with vector of pairs
  {
    Discretization d({std::make_pair(5.0, 3), std::make_pair(10.0, 5)}, 39);
    CHECK(d.getNumAdditionalPoints() == 39);
    CHECK(d.getDiscretizationPoints() == std::vector<double>{0.0, 5.0, 10.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0});
  }

  // Construct with vector of invalid pairs
  {
    CHECK_THROWS_WITH(Discretization({std::make_pair(-5.0, 3), std::make_pair(10.0, 5)}, 12),
                      Catch::Contains("Expected pairs of form [val (double), num (int)]"));
  }

  // Valid file
  {
    const Discretization d(PREPARE_DECODING_TEST_DIR "/data/test.disc");
    CHECK(d.isFile());
    CHECK(Catch::contains(d.getDiscretizationFile(), "test.disc"));
  }

  // Invalid file
  {
    CHECK_THROWS_WITH(Discretization(PREPARE_DECODING_TEST_DIR "/data/test.doesNotExist"),
                      Catch::Contains("test.doesNotExist"));
  }
}

} // namespace asmc
