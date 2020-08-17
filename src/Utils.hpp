// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#ifndef PREPAREDECODING_UTILS_HPP
#define PREPAREDECODING_UTILS_HPP

namespace asmc {

double hypergeometricPmf(int populationSize, int numberOfSuccesses, int sampleSize, int observedSuccesses);

} // namespace asmc

#endif // PREPAREDECODING_UTILS_HPP
