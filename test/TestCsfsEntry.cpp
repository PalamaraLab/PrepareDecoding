// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include "CsfsEntry.hpp"
#include "EigenTypes.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <fstream>

using Catch::Matchers::ContainsSubstring;

namespace asmc {

TEST_CASE("CSFSEntry constructor", "[CSFSEntry]") {

  std::vector<double> timeVector(5);
  std::vector<double> sizeVector(5);
  std::generate(timeVector.begin(), timeVector.end(), std::rand);
  std::generate(sizeVector.begin(), sizeVector.end(), std::rand);
  const double mu = 1.23;
  const double from = 2.34;
  const double to = 3.45;
  const int samples = 7;
  const mat_dt csfs = mat_dt::Random(3, samples - 1);

  CHECK_NOTHROW(CSFSEntry(timeVector, sizeVector, mu, from, to, samples, csfs));
}

TEST_CASE("CSFSEntry constructor throws on array size mismatch", "[CSFSEntry]") {

  std::vector<double> timeVector(5);
  std::vector<double> sizeVector(6);
  std::generate(timeVector.begin(), timeVector.end(), std::rand);
  std::generate(sizeVector.begin(), sizeVector.end(), std::rand);
  const double mu = 1.23;
  const double from = 2.34;
  const double to = 3.45;
  const int samples = 7;
  const mat_dt csfs = mat_dt::Random(3, samples - 1);

  CHECK_THROWS_WITH(CSFSEntry(timeVector, sizeVector, mu, from, to, samples, csfs),
                    ContainsSubstring("Time vector:\n"));
}

TEST_CASE("CSFSEntry constructor throws on wrong CSFS size", "[CSFSEntry]") {

  std::vector<double> timeVector(5);
  std::vector<double> sizeVector(5);
  std::generate(timeVector.begin(), timeVector.end(), std::rand);
  std::generate(sizeVector.begin(), sizeVector.end(), std::rand);
  const double mu = 1.23;
  const double from = 2.34;
  const double to = 3.45;
  const int samples = 7;

  CHECK_THROWS_WITH(CSFSEntry(timeVector, sizeVector, mu, from, to, samples, mat_dt::Random(2, samples - 1)),
                    ContainsSubstring("Time vector:\n"));

  CHECK_THROWS_WITH(CSFSEntry(timeVector, sizeVector, mu, from, to, samples, mat_dt::Random(3, samples - 2)),
                    ContainsSubstring("Time vector:\n"));
}

TEST_CASE("CSFSEntry constructor throws with from >= to", "[CSFSEntry]") {

  std::vector<double> timeVector(5);
  std::vector<double> sizeVector(5);
  std::generate(timeVector.begin(), timeVector.end(), std::rand);
  std::generate(sizeVector.begin(), sizeVector.end(), std::rand);
  const double mu = 1.23;
  const double from = 4.56;
  const double to = 3.45;
  const int samples = 7;
  const mat_dt csfs = mat_dt::Random(3, samples - 1);

  CHECK_THROWS_WITH(CSFSEntry(timeVector, sizeVector, mu, from, to, samples, csfs),
                    ContainsSubstring("Time vector:\n"));
}

TEST_CASE("CSFSEntry toString method", "[CSFSEntry]") {

  std::vector<double> timeVector(3);
  std::vector<double> sizeVector(3);
  std::fill(timeVector.begin(), timeVector.end(), 1.2);
  std::fill(sizeVector.begin(), sizeVector.end(), 2.3);
  const double mu = 3.4;
  const double from = 4.5;
  const double to = 5.6;
  const int samples = 7;
  const mat_dt csfs = mat_dt::Ones(3, samples - 1) * 8.9;

  auto csfs_entry = CSFSEntry(timeVector, sizeVector, mu, from, to, samples, csfs);

  REQUIRE(csfs_entry.toString() ==
          "Time:\t1.2 1.2 1.2\nSize:\t2.3 2.3 2.3\nMu:\t3.4\nSamples:\t7\nInterval:\t4.5\t5.6\n8.9 8.9 8.9 8.9 8.9 "
          "8.9\n8.9 8.9 8.9 8.9 8.9 8.9\n8.9 8.9 8.9 8.9 8.9 8.9");
}

} // namespace asmc
