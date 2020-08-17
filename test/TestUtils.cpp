// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include "Utils.hpp"

#include <catch2/catch.hpp>

namespace asmc {

TEST_CASE("Util: hypergeometric PMF", "[Utils]") {

  Catch::StringMaker<double>::precision = 18;

  CHECK(hypergeometricPmf(10'000, 4'270, 300, 87) == Approx(2.209102501827205049068452429E-7).epsilon(1e-12));
  CHECK(hypergeometricPmf(10'000, 4'270, 300, 128) == Approx(0.047240111763330826329736938358).epsilon(1e-12));

}

} // namespace asmc
