// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include "library.hpp"

#include <catch2/catch.hpp>

#include <filesystem>

TEST_CASE("Require true", "[test_tag]") {
  REQUIRE(true);
}

TEST_CASE("Hello", "[test_tag]") {
  REQUIRE(hello() == "Hello, world!\n");
}

TEST_CASE("Eigen", "[test_tag]") {
  REQUIRE(doSomethingWithEigen() == 3.0);
}

TEST_CASE("Zlib", "[test_tag]") {
  doSomethingWithZlib();
  std::filesystem::path path_to_test_file{"test_gz_file.gz"};
  REQUIRE(std::filesystem::is_regular_file(path_to_test_file));
  std::filesystem::remove(path_to_test_file);
}
