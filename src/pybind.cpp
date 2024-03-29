//    This file is part of ASMC, developed by Pier Francesco Palamara.
//
//    ASMC is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    ASMC is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with ASMC.  If not, see <https://www.gnu.org/licenses/>.

#include "Csfs.hpp"
#include "CsfsEntry.hpp"
#include "DecodingQuantities.hpp"
#include "EigenTypes.hpp"
#include "PrepareDecoding.hpp"
#include "ThinParameterTypes.hpp"
#include "Transition.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <vector>

namespace py = pybind11;
using namespace py::literals;
using namespace asmc;

PYBIND11_MAKE_OPAQUE(DecodingQuantities)
PYBIND11_MAKE_OPAQUE(CSFSEntry)
PYBIND11_MAKE_OPAQUE(std::vector<double>)
PYBIND11_MAKE_OPAQUE(std::vector<int>)
PYBIND11_MAKE_OPAQUE(std::vector<vec_dt>)
PYBIND11_MAKE_OPAQUE(std::vector<mat_dt>)
PYBIND11_MAKE_OPAQUE(std::map<double, CSFSEntry>)

PYBIND11_MODULE(preparedecoding_python_bindings, m) {
    py::enum_<TransitionType>(m, "TransitionType", py::arithmetic())
        .value("SMC", TransitionType::SMC)
        .value("SMC1", TransitionType::SMC1)
        .value("CSC", TransitionType::CSC);
    py::bind_vector<std::vector<int>>(m, "VectorInt");
    py::bind_vector<std::vector<double>>(m, "VectorDouble");
    py::bind_vector<std::vector<vec_dt>>(m, "VectorEigenVector");
    py::bind_vector<std::vector<mat_dt>>(m, "VectorEigenMatrix");
    py::bind_map<std::map<double, CSFSEntry>>(m, "MapDoubleToCSFSEntry");
    py::class_<DecodingQuantities>(m, "DecodingQuantities")
        .def_property_readonly("times", &DecodingQuantities::getTimes)
        .def_property_readonly("sizes", &DecodingQuantities::getSizes)
        .def("discretization", &DecodingQuantities::getDiscretization)
        .def_property_readonly("expectedTimes", &DecodingQuantities::getExpectedTimes)
        .def_property_readonly("geneticDistances", &DecodingQuantities::getGeneticDistances)
        .def_property_readonly("physDistances", &DecodingQuantities::getPhysDistances)
        .def_property_readonly("Dvectors", &DecodingQuantities::getDvectors)
        .def_property_readonly("Bvectors", &DecodingQuantities::getBvectors)
        .def_property_readonly("Uvectors", &DecodingQuantities::getBvectors)
        .def_property_readonly("rowRatioVectors", &DecodingQuantities::getRowRatioVectors)
        .def_property_readonly("homozygousEmissions", &DecodingQuantities::getHomozygousEmissions)
        .def_property_readonly("columnRatios", &DecodingQuantities::getColumnRatios)
        .def_property_readonly("states", &DecodingQuantities::getStates)
        .def_property_readonly("samples", &DecodingQuantities::getCSFSSamples)
        .def_property_readonly("mu", &DecodingQuantities::getMu)
        .def_property_readonly("CSFS", &DecodingQuantities::getCSFS)
        .def_property_readonly("foldedCSFS", &DecodingQuantities::getFoldedCSFS)
        .def_property_readonly("ascertainedCSFS", &DecodingQuantities::getAscertainedCSFS)
        .def_property_readonly("foldedAscertainedCSFS", &DecodingQuantities::getFoldedAscertainedCSFS)
        .def_property_readonly("classicEmission", &DecodingQuantities::getClassicEmission)
        .def_property_readonly("compressedEmission", &DecodingQuantities::getCompressedEmission)
        .def("save_intervals", &DecodingQuantities::saveIntervals, "Save DecodingQuantities intervals",
             "outputFileRoot"_a)
        .def("save_decoding_quantities", &DecodingQuantities::saveDecodingQuantities, "Save DecodingQuantities",
             "outputFileRoot"_a)
        .def("save_csfs", &DecodingQuantities::saveCsfs, "Save CSFS",
             "outputFileRoot"_a)
        .def("save_discretization", &DecodingQuantities::saveDiscretization, "Save discretization",
             "outputFileRoot"_a);
    py::class_<CSFS>(m, "CSFS")
      .def_static("load", &CSFS::load, "Construct CSFS from parameters",
          "times"_a, "sizes"_a, "mu"_a, "samples"_a,
          "froms"_a, "tos"_a, "csfses"_a)
      .def_static("loadFromFile", &CSFS::loadFromFile, "Load precomputed CSFS from file",
          "filename"_a);
    py::class_<CSFSEntry>(m, "CSFSEntry")
      .def_property_readonly("times", &CSFSEntry::getTime)
      .def_property_readonly("sizes", &CSFSEntry::getSize)
      .def_property_readonly("mu", &CSFSEntry::getMu)
      .def_property_readonly("samples", &CSFSEntry::getSamples)
      .def_property_readonly("csfsFrom", &CSFSEntry::getFrom)
      .def_property_readonly("csfsTo", &CSFSEntry::getTo)
      .def_property_readonly("csfs", &CSFSEntry::getCSFSMatrix)
    ;
    py::class_<Demography>(m, "Demography")
        .def(py::init<>())
        .def(py::init<std::string_view>());

    py::class_<Discretization>(m, "Discretization")
        .def(py::init<std::string_view>())
        .def(py::init<std::vector<double>, int>())
        .def(py::init<const std::vector<std::pair<double, int>>&, int>());

    py::class_<Frequencies>(m, "Frequencies")
        .def(py::init<std::string_view>())
        .def(py::init<std::string_view, unsigned>());

    m.def("prepareDecoding", &prepareDecoding, "Calculate decoding quantities", "demography"_a, "discretization"_a,
          "frequencies"_a, "csfs_file"_a = std::string_view{}, "file_root"_a = std::string_view{}, "mut_rate"_a = 1.65e-8,
          "samples"_a = 300);

    m.def("save_demography", &demo::saveDemography, "Save a default demography to file", "output_dir"_a, "demo_code"_a);
}
