// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include "CsfsEntry.hpp"
#include "EigenTypes.hpp"

#include <catch2/catch.hpp>

#include <fstream>

namespace asmc {

TEST_CASE("CsfsEntry constructor", "[CsfsEntry]") {

  const array_dt timeVector = array_dt::Random(5);
  const array_dt sizeVector = array_dt::Random(5);
  const double mu = 1.23;
  const double from = 2.34;
  const double to = 3.45;
  const int samples = 7;
  const mat_dt csfs = mat_dt::Random(3, samples - 1);

  CHECK_NOTHROW(CsfsEntry(timeVector, sizeVector, mu, from, to, samples, csfs));
}

TEST_CASE("CsfsEntry constructor throws on array size mismatch", "[CsfsEntry]") {

  const array_dt timeVector = array_dt::Random(5);
  const array_dt sizeVector = array_dt::Random(6);
  const double mu = 1.23;
  const double from = 2.34;
  const double to = 3.45;
  const int samples = 7;
  const mat_dt csfs = mat_dt::Random(3, samples - 1);

  CHECK_THROWS_WITH(CsfsEntry(timeVector, sizeVector, mu, from, to, samples, csfs),
                    Catch::Contains("Time vector:\n"));
}

TEST_CASE("CsfsEntry constructor throws on wrong CSFS size", "[CsfsEntry]") {

  const array_dt timeVector = array_dt::Random(5);
  const array_dt sizeVector = array_dt::Random(5);
  const double mu = 1.23;
  const double from = 2.34;
  const double to = 3.45;
  const int samples = 7;

  CHECK_THROWS_WITH(CsfsEntry(timeVector, sizeVector, mu, from, to, samples, mat_dt::Random(2, samples - 1)),
                    Catch::Contains("Time vector:\n"));

  CHECK_THROWS_WITH(CsfsEntry(timeVector, sizeVector, mu, from, to, samples, mat_dt::Random(3, samples - 2)),
                    Catch::Contains("Time vector:\n"));
}

TEST_CASE("CsfsEntry constructor throws with to >= from", "[CsfsEntry]") {

  const array_dt timeVector = array_dt::Random(5);
  const array_dt sizeVector = array_dt::Random(5);
  const double mu = 1.23;
  const double from = 4.56;
  const double to = 3.45;
  const int samples = 7;
  const mat_dt csfs = mat_dt::Random(3, samples - 1);

  CHECK_THROWS_WITH(CsfsEntry(timeVector, sizeVector, mu, from, to, samples, csfs),
                    Catch::Contains("Time vector:\n"));
}

TEST_CASE("CsfsEntry toString method", "[CsfsEntry]") {

  const array_dt timeVector = array_dt::Zero(3) + 1.2;
  const array_dt sizeVector = array_dt::Zero(3) + 2.3;
  const double mu = 3.4;
  const double from = 4.5;
  const double to = 5.6;
  const int samples = 7;
  const mat_dt csfs = mat_dt::Ones(3, samples - 1) * 8.9;

  auto csfs_entry = CsfsEntry(timeVector, sizeVector, mu, from, to, samples, csfs);

  REQUIRE(csfs_entry.toString() ==
          "Time:\t1.2 1.2 1.2\nSize:\t2.3 2.3 2.3\nMu:\t3.4\nSamples:\t7\nInterval:\t4.5\t5.6\n8.9 8.9 8.9 8.9 8.9 "
          "8.9\n8.9 8.9 8.9 8.9 8.9 8.9\n8.9 8.9 8.9 8.9 8.9 8.9\n");
}

} // namespace asmc