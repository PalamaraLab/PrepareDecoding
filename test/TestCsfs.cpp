// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include "Csfs.hpp"

#include <catch2/catch.hpp>

#include <fstream>

namespace asmc {

TEST_CASE("Csfs constructor", "[Csfs]") {
  Csfs csfs;
  csfs = Csfs();
}

} // namespace asmc
