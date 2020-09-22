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
    std::vector<double> mSpectrum;
    // probability of monomorphic is stored separately
    double mMonomorphicProbability;
    void normalize(std::vector<double>& spectrum);

  public:
    explicit ArraySpectrum(Data data, unsigned samples);
    double getMonomorphic() const { return mMonomorphicProbability; }

};

} // namespace asmc

#endif // PREPAREDECODING_ARRAYSPECTRUM_HPP
