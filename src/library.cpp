// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include "library.hpp"

#include <Eigen/Dense>
#include <fmt/core.h>
#include <zlib.h>

std::string hello() {
  return fmt::format("Hello, world!\n");
}

double doSomethingWithEigen() {

  using mat_dt = Eigen::MatrixXd;

  mat_dt m(2, 2);
  m(0, 0) = 3.0;

  return m(0, 0);
}

void doSomethingWithZlib() {

  auto gz_file = gzopen("test_gz_file.gz", "w");

  std::string payload = fmt::format("{}\t{}\t{}\t{}\t{}", 123, 234, 345, 456, 567);
  gzwrite(gz_file, payload.c_str(), static_cast<unsigned>(payload.size()));

  gzclose(gz_file);
}
