#define FMT_HEADER_ONLY
#include "observable.h"
#include "propagator.h"
#include "utility/abort.h"
#include "utility/fmt/format.h"
#include "utility/fmt/printf.h"
#include "utility/utility.h"
#include <cmath>
#include <iostream>

using namespace std;
using namespace obs;

extern parameter Para;
extern variable Var;

oneBodyObs::oneBodyObs() {
  Normalization = 1.0e-10;

  // The zeroth order of polar, sigma and delta all have one external K and one
  // external Tau
  if (DiagType == POLAR) {
    Name = "polar";
    ksize = Para.BoseKGrid.size;
  } else if (DiagType == SIGMA) {
    Name = "sigma";
    ksize = Para.FermiKGrid.size;
  } else if (DiagType == DELTA) {
    Name = "delta";
    ksize = Para.FermiKGrid.size;
  } else
    return;
  // ABORT("not implemented!");
  // do nothing

  PhyWeight = ksize * Para.TauGrid.size;
  _Estimator.Initialize({Para.Order + 1, ksize, Para.TauGrid.size});

  return;
}

void oneBodyObs::Measure0(double Factor) { Normalization += 1.0 * Factor; }
void oneBodyObs::Measure(int Order, int KBin, int TauBin, double Weight,
                         double Factor) {
  ASSERT(KBin >= 0 && KBin < ksize, "Kidx is out of range!");
  ASSERT(TauBin >= 0 && TauBin < Para.TauGrid.size, "TauIdx is out of range!");

  _Estimator(Order, KBin, TauBin) += Weight * Factor;
  _Estimator(0, KBin, TauBin) += Weight * Factor;
}

void oneBodyObs::Save() {

  string FileName = fmt::format("{0}_pid{1}.dat", Name, Para.PID);
  ofstream VerFile;
  VerFile.open(FileName, ios::out | ios::trunc);

  if (VerFile.is_open()) {

    VerFile << "# Counter: " << Var.Counter << endl;
    VerFile << "# Norm: " << Normalization << endl;

    if (DiagType == POLAR)
      VerFile << "# KGrid: " << Para.BoseKGrid.str() << endl;
    else
      VerFile << "# KGrid: " << Para.FermiKGrid.str() << endl;

    VerFile << "# TauGrid: " << Para.TauGrid.str() << endl;
    for (int order = 0; order <= Para.Order; order++)
      for (int qindex = 0; qindex < ksize; ++qindex)
        for (int tindex = 0; tindex < Para.TauGrid.size; ++tindex)
          VerFile << _Estimator(order, qindex, tindex) * PhyWeight << "  ";
    VerFile.close();
  } else {
    LOG_WARNING(Name << " for PID " << Para.PID << " fails to save!");
  }
}

ver4Obs::ver4Obs() {
  Normalization = 1.0e-10;
  PhyWeight = Para.AngleGrid.size * Para.Beta;
  for (auto &estimator : _Estimator)
    estimator.Initialize(
        {Para.Order + 1, Para.AngleGrid.size, Para.BoseKGrid.size});
};

void ver4Obs::Measure0(double Factor) { Normalization += 1.0 * Factor; }

void ver4Obs::Measure(int Order, int QIndex, int AngleIndex,
                      const std::vector<verWeight> &Weight, double Factor) {

  ASSERT(Order != 0, "Order must be >=1!");

  ASSERT(AngleIndex >= 0 && AngleIndex < Para.AngleGrid.size,
         "AngleIndex out of range!");

  for (int chan = 0; chan < 4; ++chan) {
    // cout << "chan=" << chan << ", " << AngleIndex << ", " << QIndex << endl;
    // cout << Weight[chan] << endl;
    _Estimator[chan](Order, AngleIndex, QIndex) += Weight[chan] * Factor;
    _Estimator[chan](0, AngleIndex, QIndex) += Weight[chan] * Factor;
  }
  return;
}
void ver4Obs::Save() {
  string FileName = fmt::format("vertex_pid{0}.dat", Para.PID);
  ofstream VerFile;
  VerFile.open(FileName, ios::out | ios::trunc);

  if (VerFile.is_open()) {

    VerFile << "# Counter: " << Var.Counter << endl;
    VerFile << "# Norm: " << Normalization << endl;
    VerFile << "# KGrid: " << Para.BoseKGrid.str() << endl;
    VerFile << "# AngleGrid: " << Para.AngleGrid.str() << endl;

    // cout << "Size: " << Para.Order << ", " << Para.AngleGrid.size << ", "
    //      << Para.BoseKGrid.size << endl;

    for (int order = 0; order <= Para.Order; order++)
      for (int chan = 0; chan < 4; chan++)
        for (int angle = 0; angle < Para.AngleGrid.size; ++angle)
          for (int qindex = 0; qindex < Para.BoseKGrid.size; ++qindex)
            for (int dir = 0; dir < 2; ++dir) {
              VerFile << _Estimator[chan](order, angle, qindex)[dir] * PhyWeight
                      << "  ";
            }
    VerFile.close();
  } else
    LOG_WARNING("Vertex4 for PID " << Para.PID << " fails to save!");
}