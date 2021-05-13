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
    const Demography d;
    CHECK_THROWS_WITH(Demography(PREPARE_DECODING_TEST_DIR "/data/test.doesNotExist"),
                      Catch::Contains("test.doesNotExist"));
  }
}

} // namespace asmc
