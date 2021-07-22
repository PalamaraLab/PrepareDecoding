// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include "Data.hpp"
#include "ThinParameterTypes.hpp"

#include <catch2/catch.hpp>

namespace asmc {

TEST_CASE("Data explicit constructor, no hap file found", "[Data]") {

  CHECK_THROWS_WITH(Data("some/path"), Catch::Contains("No haps file found at some/path"));

//  Data data("test/data/data");
}

TEST_CASE("Data read frequency file", "[Data]") {

   std::vector<double> freq{0.1684, 0.03716, 0.1382, 0.1465, 0.1098};
   std::vector<unsigned int> counts{302964, 304088, 303726, 287624, 303648};
   std::vector<unsigned int> minorAlleles{51019, 11299, 41974, 42136, 33340};

   {
     Data data;
     data.addFreq(Frequencies(PREPARE_DECODING_TEST_DIR "/data/example.frq"));
     CHECK(data.getAllSNPsFreq() == freq);
     CHECK(data.getAllSNPsAlleleCounts() == counts);
     CHECK(data.getAllSNPsMinorAlleles() == minorAlleles);
   }

}

} // namespace asmc
