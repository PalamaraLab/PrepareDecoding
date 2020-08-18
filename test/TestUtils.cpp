// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include "Utils.hpp"

#include <catch2/catch.hpp>

namespace asmc {

TEST_CASE("Util: hypergeometric PMF specific values", "[Utils]") {

  Catch::StringMaker<double>::precision = 18;

  CHECK(hypergeometricPmf(10'000, 4'270, 300, 87) == Approx(2.209102501827205049068452429E-7).epsilon(1e-9));
  CHECK(hypergeometricPmf(10'000, 4'270, 300, 128) == Approx(0.047240111763330826329736938358).epsilon(1e-9));
}

TEST_CASE("Util: hypergeometric PMF sums to one", "[Utils]") {

  Catch::StringMaker<double>::precision = 18;

  const int populationSize = 12'345;
  const int numberOfSuccesses = 6'789;
  const int sampleSize = 300;

  double total_probability = 0.0;
  for (int observedSuccesses = 0; observedSuccesses < sampleSize; ++observedSuccesses) {
    total_probability += hypergeometricPmf(populationSize, numberOfSuccesses, sampleSize, observedSuccesses);
  }

  CHECK(total_probability == Approx(1.0).epsilon(1e-9));
}

} // namespace asmc
