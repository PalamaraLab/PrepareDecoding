// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include "library.hpp"

#include <catch2/catch.hpp>

TEST_CASE("Require true", "[test_tag]") {
  REQUIRE(true);
}

TEST_CASE("Hello", "[test_tag]") {
  REQUIRE(hello() == "Hello, world!\n");
}

TEST_CASE("Eigen", "[test_tag]") {
  REQUIRE(doSomethingWithEigen() == 3.0);
}
