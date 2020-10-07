#include "../src/lib/grid.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace grid;

PYBIND11_MODULE(grid, m) {
  py::class_<Tau>(m, "Tau")
      // .def("Size", &tauGrid::Size)
      .def(py::init<>())
      .def("build", &Tau::build)
      .def("floor", &Tau::floor)
      .def("str", &Tau::str)
      .def_readwrite("size", &Tau::size)
      .def_readwrite("grid", &Tau::grid)
      .def_readwrite("weight", &Tau::weight);

  py::class_<TauUL>(m, "TauUL")
    // .def("Size", &tauGrid::Size)
    .def(py::init<>())
    .def("build", &TauUL::build)
    .def("floor", &TauUL::floor)
    .def("str", &TauUL::str)
    .def_readwrite("size", &TauUL::size)
    .def_readwrite("grid", &TauUL::grid)
    .def_readwrite("weight", &TauUL::weight);

  py::class_<FermiK>(m, "FermiK")
      // .def("Size", &tauGrid::Size)
      .def(py::init<>())
      .def("build", &FermiK::build)
      .def("floor", &FermiK::floor)
      .def("str", &FermiK::str)
      .def_readwrite("size", &FermiK::size)
      .def_readwrite("kFidx", &FermiK::kFidx)
      .def_readwrite("grid", &FermiK::grid);

  py::class_<FermiKUL>(m, "FermiKUL")
      // .def("Size", &tauGrid::Size)
      .def(py::init<>())
      .def("build", &FermiKUL::build)
      .def("floor", &FermiKUL::floor)
      .def("str", &FermiKUL::str)
      .def_readwrite("size", &FermiKUL::size)
      .def_readwrite("kFidx", &FermiKUL::kFidx)
      .def_readwrite("grid", &FermiKUL::grid);

  py::class_<BoseK>(m, "BoseK")
      // .def("Size", &tauGrid::Size)
      .def(py::init<>())
      .def("build", &BoseK::build)
      .def("floor", &BoseK::floor)
      .def("str", &BoseK::str)
      .def_readwrite("size", &BoseK::size)
      .def_readwrite("kFidx", &BoseK::kFidx)
      .def_readwrite("twokFidx", &BoseK::twokFidx)
      .def_readwrite("grid", &BoseK::grid);

  py::class_<BoseKUL>(m, "BoseKUL")
      // .def("Size", &tauGrid::Size)
      .def(py::init<>())
      .def("build", &BoseKUL::build)
      .def("floor", &BoseKUL::floor)
      .def("str", &BoseKUL::str)
      .def_readwrite("size", &BoseKUL::size)
      .def_readwrite("kFidx", &BoseKUL::kFidx)
      .def_readwrite("twokFidx", &BoseKUL::twokFidx)
      .def_readwrite("grid", &BoseKUL::grid);

  py::class_<Uniform>(m, "Uniform")
      // .def("Size", &tauGrid::Size)
      .def(py::init<>())
      .def("build", &Uniform::build)
      .def("floor", &Uniform::floor)
      .def("str", &Uniform::str)
      .def_readwrite("size", &Uniform::size)
      .def_readwrite("grid", &Uniform::grid);
}
