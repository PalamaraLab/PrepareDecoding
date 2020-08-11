// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include "library.hpp"

#include <Eigen/Dense>
#include <fmt/core.h>

std::string hello() {
  return fmt::format("Hello, world!\n");
}

double doSomethingWithEigen() {

  using mat_t = Eigen::MatrixXd;

  mat_t m(2, 2);
  m(0, 0) = 3.0;

  return m(0, 0);
}
