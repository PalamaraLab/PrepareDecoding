// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include "CsfsEntry.hpp"

#include <catch2/catch.hpp>

using array_t = Eigen::ArrayXd;
using mat_t = Eigen::MatrixXd;

TEST_CASE("Can successfully construct object", "[CsfsEntry]") {

  const array_t timeVector = array_t::Random(5);
  const array_t sizeVector = array_t::Random(5);
  const double mu = 1.23;
  const double from = 2.34;
  const double to = 3.45;
  const int samples = 7;
  const mat_t csfs = mat_t::Random(3, samples - 1);

  CHECK_NOTHROW(asmc::CsfsEntry(timeVector, sizeVector, mu, from, to, samples, csfs));
}

TEST_CASE("Throws on array size mismatch", "[CsfsEntry]") {

  const array_t timeVector = array_t::Random(5);
  const array_t sizeVector = array_t::Random(6);
  const double mu = 1.23;
  const double from = 2.34;
  const double to = 3.45;
  const int samples = 7;
  const mat_t csfs = mat_t::Random(3, samples - 1);

  CHECK_THROWS_WITH(asmc::CsfsEntry(timeVector, sizeVector, mu, from, to, samples, csfs),
                    Catch::Contains("Time vector:\n"));
}

TEST_CASE("Throws on wrong CSFS size", "[CsfsEntry]") {

  const array_t timeVector = array_t::Random(5);
  const array_t sizeVector = array_t::Random(5);
  const double mu = 1.23;
  const double from = 2.34;
  const double to = 3.45;
  const int samples = 7;

  CHECK_THROWS_WITH(asmc::CsfsEntry(timeVector, sizeVector, mu, from, to, samples, mat_t::Random(2, samples - 1)),
                    Catch::Contains("Time vector:\n"));

  CHECK_THROWS_WITH(asmc::CsfsEntry(timeVector, sizeVector, mu, from, to, samples, mat_t::Random(3, samples - 2)),
                    Catch::Contains("Time vector:\n"));
}

TEST_CASE("Throws with to >= from", "[CsfsEntry]") {

  const array_t timeVector = array_t::Random(5);
  const array_t sizeVector = array_t::Random(5);
  const double mu = 1.23;
  const double from = 4.56;
  const double to = 3.45;
  const int samples = 7;
  const mat_t csfs = mat_t::Random(3, samples - 1);

  CHECK_THROWS_WITH(asmc::CsfsEntry(timeVector, sizeVector, mu, from, to, samples, csfs),
                    Catch::Contains("Time vector:\n"));
}
