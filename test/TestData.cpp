// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include "Data.hpp"

#include <catch2/catch.hpp>

namespace asmc {

TEST_CASE("Data default constructor", "[Data]") {
  Data data;
}

TEST_CASE("Data explicit constructor, no hap file found", "[Data]") {

  CHECK_THROWS_WITH(Data("some/path"), Catch::Contains("No haps file found at some/path"));

//  Data data("test/data/data");
}

} // namespace asmc
