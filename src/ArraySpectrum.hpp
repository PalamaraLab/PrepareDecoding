// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include <vector>
#include "Data.hpp"

#ifndef PREPAREDECODING_ARRAYSPECTRUM_HPP
#define PREPAREDECODING_ARRAYSPECTRUM_HPP

namespace asmc {

class ArraySpectrum {

  private:
    // spectrum does not include probability of monomorphic alleles due
    // to no variation or subsampling
    array_dt mSpectrum = {};
    // probability of monomorphic is stored separately
    double mMonomorphicProbability = 0.0;

  public:
    ArraySpectrum() = default;
    explicit ArraySpectrum(Data data, unsigned samples);
    array_dt getSpectrum() { return mSpectrum; }
    double getMonomorphic() const { return mMonomorphicProbability; }

};

} // namespace asmc

#endif // PREPAREDECODING_ARRAYSPECTRUM_HPP
