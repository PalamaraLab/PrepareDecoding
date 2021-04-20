// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#ifndef PREPAREDECODING_SMCPP_HPP
#define PREPAREDECODING_SMCPP_HPP

#include <Eigen/Dense>
#include <vector>

template <typename T> using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template <typename T> using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

/**
 * Initialise cache to store smcpp binary data objects. This will be in a temporary directory.
 */
void smcpp_init_cache();

/**
 * Compute the raw sfs, and return it as a Matrix of doubles
 *
 * @param a normalised population sizes from the demographic history
 * @param s discrete derivative of the times at which demographic history population sizes were sampled
 * @param n number of samples
 * @param t1 start time (in generations) of the discrete time interval
 * @param t2 end time (in generations) of the discrete time interval
 * @param below_only whether to only compute below
 * @return the calculated raw sfs
 */
Matrix<double> raw_sfs(const std::vector<double>& a, const std::vector<double>& s, int n, double t1, double t2,
                       bool below_only = false);

#endif // PREPAREDECODING_SMCPP_HPP
