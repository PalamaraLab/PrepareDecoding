// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#ifndef PREPAREDECODING_UTILS_HPP
#define PREPAREDECODING_UTILS_HPP

#include <zlib.h>
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>

namespace asmc {

double hypergeometricPmf(int populationSize, int numberOfSuccesses, int sampleSize, int observedSuccesses);

/**
 * Read the next line from a gzip file.
 *
 * @param gzFileHandle handle to a file opened with zlib's gzopen
 * @return a string contating the next line contained in the gzip file, without a trailing newline character
 */
std::string readNextLineFromGzip(gzFile& gzFileHandle);

void normalize(std::vector<double>& spectrum);

} // namespace asmc

#endif // PREPAREDECODING_UTILS_HPP
