#include "../src/lib/green.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
PYBIND11_MODULE(green, m) {
  m.def("fermiGreen", &fermiGreen, "bare fermionic Green's function");
  m.def("fockYukawa", &fockYukawa,
        "zero temperature Yukawa interaction Fock diagram");
}
