// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include <cxxopts.hpp>
#include <fmt/core.h>

int main(int argc, char* argv[]) {

  cxxopts::Options options("MyProgram", "One line description of MyProgram");

  options.add_options()
          ("i,integer", "Int param", cxxopts::value<int>()->default_value("-1"))
          ("f,file", "File name", cxxopts::value<std::string>()->default_value("/some/file"))
          ("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"));

  auto result = options.parse(argc, argv);

  fmt::print("Option i: {}\n", result["i"].as<int>());
  fmt::print("Option f: {}\n", result["f"].as<std::string>());
  fmt::print("Option v: {}\n", result["v"].as<bool>());
}