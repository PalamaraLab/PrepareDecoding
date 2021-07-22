// This file is part of https://github.com/PalamaraLab/PrepareDecoding which is released under the GPL-3.0 license.
// See accompanying LICENSE and COPYING for copyright notice and full details.

#include "DefaultDemographies.hpp"

#include <fmt/format.h>
#include <fmt/os.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <filesystem>
#include <string_view>
#include <tuple>
#include <vector>

namespace fs = std::filesystem;

namespace asmc::demo {

std::tuple<std::vector<double>, std::vector<double>> getBuiltInDemography(std::string_view demoCode) {
  std::vector<double> times;
  std::vector<double> sizes;
  if (demoCode == "ACB") {
    times = std::vector<double>(timesACB.begin(), timesACB.end());
    sizes = std::vector<double>(sizesACB.begin(), sizesACB.end());
  } else if (demoCode == "ASW") {
    times = std::vector<double>(timesASW.begin(), timesASW.end());
    sizes = std::vector<double>(sizesASW.begin(), sizesASW.end());
  } else if (demoCode == "BEB") {
    times = std::vector<double>(timesBEB.begin(), timesBEB.end());
    sizes = std::vector<double>(sizesBEB.begin(), sizesBEB.end());
  } else if (demoCode == "CDX") {
    times = std::vector<double>(timesCDX.begin(), timesCDX.end());
    sizes = std::vector<double>(sizesCDX.begin(), sizesCDX.end());
  } else if (demoCode == "CEU") {
    times = std::vector<double>(timesCEU.begin(), timesCEU.end());
    sizes = std::vector<double>(sizesCEU.begin(), sizesCEU.end());
  } else if (demoCode == "CHB") {
    times = std::vector<double>(timesCHB.begin(), timesCHB.end());
    sizes = std::vector<double>(sizesCHB.begin(), sizesCHB.end());
  } else if (demoCode == "CHS") {
    times = std::vector<double>(timesCHS.begin(), timesCHS.end());
    sizes = std::vector<double>(sizesCHS.begin(), sizesCHS.end());
  } else if (demoCode == "CLM") {
    times = std::vector<double>(timesCLM.begin(), timesCLM.end());
    sizes = std::vector<double>(sizesCLM.begin(), sizesCLM.end());
  } else if (demoCode == "ESN") {
    times = std::vector<double>(timesESN.begin(), timesESN.end());
    sizes = std::vector<double>(sizesESN.begin(), sizesESN.end());
  } else if (demoCode == "FIN") {
    times = std::vector<double>(timesFIN.begin(), timesFIN.end());
    sizes = std::vector<double>(sizesFIN.begin(), sizesFIN.end());
  } else if (demoCode == "GBR") {
    times = std::vector<double>(timesGBR.begin(), timesGBR.end());
    sizes = std::vector<double>(sizesGBR.begin(), sizesGBR.end());
  } else if (demoCode == "GIH") {
    times = std::vector<double>(timesGIH.begin(), timesGIH.end());
    sizes = std::vector<double>(sizesGIH.begin(), sizesGIH.end());
  } else if (demoCode == "GWD") {
    times = std::vector<double>(timesGWD.begin(), timesGWD.end());
    sizes = std::vector<double>(sizesGWD.begin(), sizesGWD.end());
  } else if (demoCode == "IBS") {
    times = std::vector<double>(timesIBS.begin(), timesIBS.end());
    sizes = std::vector<double>(sizesIBS.begin(), sizesIBS.end());
  } else if (demoCode == "ITU") {
    times = std::vector<double>(timesITU.begin(), timesITU.end());
    sizes = std::vector<double>(sizesITU.begin(), sizesITU.end());
  } else if (demoCode == "JPT") {
    times = std::vector<double>(timesJPT.begin(), timesJPT.end());
    sizes = std::vector<double>(sizesJPT.begin(), sizesJPT.end());
  } else if (demoCode == "KHV") {
    times = std::vector<double>(timesKHV.begin(), timesKHV.end());
    sizes = std::vector<double>(sizesKHV.begin(), sizesKHV.end());
  } else if (demoCode == "LWK") {
    times = std::vector<double>(timesLWK.begin(), timesLWK.end());
    sizes = std::vector<double>(sizesLWK.begin(), sizesLWK.end());
  } else if (demoCode == "MSL") {
    times = std::vector<double>(timesMSL.begin(), timesMSL.end());
    sizes = std::vector<double>(sizesMSL.begin(), sizesMSL.end());
  } else if (demoCode == "MXL") {
    times = std::vector<double>(timesMXL.begin(), timesMXL.end());
    sizes = std::vector<double>(sizesMXL.begin(), sizesMXL.end());
  } else if (demoCode == "PEL") {
    times = std::vector<double>(timesPEL.begin(), timesPEL.end());
    sizes = std::vector<double>(sizesPEL.begin(), sizesPEL.end());
  } else if (demoCode == "PJL") {
    times = std::vector<double>(timesPJL.begin(), timesPJL.end());
    sizes = std::vector<double>(sizesPJL.begin(), sizesPJL.end());
  } else if (demoCode == "PUR") {
    times = std::vector<double>(timesPUR.begin(), timesPUR.end());
    sizes = std::vector<double>(sizesPUR.begin(), sizesPUR.end());
  } else if (demoCode == "STU") {
    times = std::vector<double>(timesSTU.begin(), timesSTU.end());
    sizes = std::vector<double>(sizesSTU.begin(), sizesSTU.end());
  } else if (demoCode == "TSI") {
    times = std::vector<double>(timesTSI.begin(), timesTSI.end());
    sizes = std::vector<double>(sizesTSI.begin(), sizesTSI.end());
  } else if (demoCode == "YRI") {
    times = std::vector<double>(timesYRI.begin(), timesYRI.end());
    sizes = std::vector<double>(sizesYRI.begin(), sizesYRI.end());
  }
  return std::make_tuple(times, sizes);
}

bool isValidDemography(std::string_view demoCode) {
  if (std::find(std::begin(validDemographies), std::end(validDemographies), demoCode) != std::end(validDemographies)) {
    return true;
  }
  return false;
}

void saveDemography(std::string_view outputDir, std::string_view demoCode) {

  if (!fs::is_directory(outputDir)) {
    throw std::runtime_error(
        fmt::format("Error saving demography {}: {} is not a valid directory\n", demoCode, outputDir));
  }

  if (!isValidDemography(demoCode)) {
    throw std::runtime_error(fmt::format("Error saving demography: {} is not a valid demography (expected one of {})\n",
                                         demoCode, validDemographies));
  }

  auto [times, sizes] = getBuiltInDemography(demoCode);
  std::string outputFilename = (fs::path(outputDir) / demoCode).string();
  auto fmtOutFile = fmt::output_file(fmt::format("{}.demo", outputFilename));

  for (auto i = 0ul; i < times.size(); ++i) {
    fmtOutFile.print("{:#}\t{:#}\n", times.at(i), sizes.at(i));
  }
}

} // namespace asmc::demo
