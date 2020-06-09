#include "./src/grid.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>

namespace py = pybind11;

PYBIND11_MODULE(grid, m) {
  py::class_<logGrid>(m, "logGrid")
      .def("Bound", &logGrid::Bound)
      .def("Idx", &logGrid::Idx)
      .def("Lambda", &logGrid::Lambda)
      .def("Initialize", &logGrid::Initialize)
      .def("Floor", &logGrid::Floor)
      .def("Grid", &logGrid::Grid);

  // py::class_<tauGrid>(m, "tauGrid")
  //     .def("Size", &tauGrid::Size)
  //     .def("Initialize", &tauGrid::Initialize)
  //     .def("Floor", &tauGrid::Floor)
  //     .def("Grid", &tauGrid::Grid)
  //     .def("Weight", &tauGrid::Weight)
  //     .def("ToString", &tauGrid::ToString);
}

// PYBIND11_MODULE(grid, m) {
//   py::class_<tauGrid>(m, "tauGrid")
//       .def("Size", &tauGrid::Size)
//       .def("Initialize", &tauGrid::Initialize)
//       .def("Floor", &tauGrid::Floor)
//       .def("Grid", &tauGrid::Grid)
//       .def("Weight", &tauGrid::Weight)
//       .def("ToString", &tauGrid::ToString);
// }
// int add(int i, int j) { return i + j; }

// std::vector<double> TGrid(double Beta, int Size, double scale) {
//   std::vector<double> Grid;
//   Grid.resize(Size);
//   double lambda = Beta / scale / (Size / 2.0);

//   logGrid _Grid0;
//   logGrid _Grid1;

//   _Grid0.Initialize({0.0, Beta / 2.0}, {0.0, Size / 2 - 0.5}, lambda, true);
//   _Grid1.Initialize({Beta / 2.0, Beta}, {Size / 2 - 0.5, Size - 1.0}, lambda,
//                     false);

//   for (int i = 0; i < Size / 2; ++i) {
//     Grid[i] = _Grid0.Grid(i);
//   }

//   for (int i = Size / 2; i < Size; ++i) {
//     Grid[i] = _Grid1.Grid(i);
//   }

//   Grid[0] = 1.0e-8;
//   Grid[Size - 1] = Beta - 1.0e-8;
//   return Grid;
// }

// PYBIND11_MODULE(grid, m) {
// m.doc() = "pybind11 example plugin"; // optional module docstring
// m.def("add", &add, "add fnction");

// m.def("TGrid", &TGrid, "taugrid function");
// }