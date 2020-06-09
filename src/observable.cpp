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
  PhyWeight = Para.KGrid.Size * Para.TauGrid.Size;

  if (DiagType == POLAR) {
    Name = "polar";
  } else if (DiagType == SIGMA) {
    Name = "sigma";
  } else if (DiagType == DELTA) {
    Name = "delta";
  } else
    // do nothing
    return;

  _Estimator.Initialize({Para.Order + 1, Para.KGrid.Size, Para.TauGrid.Size});
}

void oneBodyObs::Measure0(double Factor) { Normalization += 1.0 * Factor; }
void oneBodyObs::Measure(int Order, int KBin, int TauBin, double Weight,
                         double Factor) {
  ASSERT(KBin >= 0 && KBin < Para.KGrid.Size, "Kidx is out of range!");
  ASSERT(TauBin >= 0 && TauBin < Para.TauGrid.Size, "TauIdx is out of range!");

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
    VerFile << "# KGrid: " << Para.KGrid.ToString() << endl;
    VerFile << "# TauGrid: " << Para.TauGrid.ToString() << endl;
    for (int order = 0; order <= Para.Order; order++)
      for (int qindex = 0; qindex < Para.KGrid.Size; ++qindex)
        for (int tindex = 0; tindex < Para.TauGrid.Size; ++tindex)
          VerFile << _Estimator(order, qindex, tindex) * PhyWeight << "  ";
    VerFile.close();
  } else {
    LOG_WARNING(Name << " for PID " << Para.PID << " fails to save!");
  }
  FileName = fmt::format("{0}_pid{1}_verb.dat", Name, Para.PID);
  VerFile.open(FileName, ios::out | ios::trunc);

  if (VerFile.is_open()) {

    VerFile << "# Counter: " << Var.Counter << endl;
    VerFile << "# Norm: " << Normalization << endl;
    for (int order = 0; order <= Para.Order; order++)
      for (int qindex = 0; qindex < Para.KGrid.Size; ++qindex)
        for (int tindex = 0; tindex < Para.TauGrid.Size; ++tindex)
          VerFile << order << "\t"
                  << Para.KGrid.Grid[qindex] << "\t"
                  << Para.TauGrid.Grid[tindex] << "\t"
                  << _Estimator(order, qindex, tindex) * PhyWeight/Normalization << "\n";
    VerFile.close();
  } else {
    LOG_WARNING(Name << " for PID " << Para.PID << " fails to save!");
  }
}

void oneBodyObs::Save(int channel) {

  string FileName = fmt::format("{0}_chan{1}_pid{2}.dat", Name,channel, Para.PID);
  ofstream VerFile;
  VerFile.open(FileName, ios::out | ios::trunc);

  if (VerFile.is_open()) {

    VerFile << "# Counter: " << Var.Counter << endl;
    VerFile << "# Norm: " << Normalization << endl;
    for (int order = 0; order <= Para.Order; order++)
      for (int qindex = 0; qindex < Para.ExtMomBinSize; ++qindex)
        for (int tindex = 0; tindex < Para.TauBinSize; ++tindex)
          VerFile << _Estimator(order, qindex, tindex) * PhyWeight << "  ";
    VerFile.close();
  } else {
    LOG_WARNING(Name << " for PID " << Para.PID << " fails to save!");
  }
  FileName = fmt::format("{0}_chan{1}_pid{2}_verb.dat", Name,channel, Para.PID);
  VerFile.open(FileName, ios::out | ios::trunc);

  if (VerFile.is_open()) {

    VerFile << "# Counter: " << Var.Counter << endl;
    VerFile << "# Norm: " << Normalization << endl;
    for (int order = 0; order <= Para.Order; order++)
      for (int qindex = 0; qindex < Para.ExtMomBinSize; ++qindex)
        for (int tindex = 0; tindex < Para.TauBinSize; ++tindex)
          VerFile << order << "\t"
                  << Para.ExtMomTable[qindex].norm() << "\t"
                  << Para.ExtTauTable[tindex] << "\t"
                  << _Estimator(order, qindex, tindex) * PhyWeight/Normalization << "\n";
    VerFile.close();
  } else {
    LOG_WARNING(Name << " for PID " << Para.PID << " fails to save!");
  }
}

ver4Obs::ver4Obs() {
  Normalization = 1.0e-10;
  PhyWeight = Para.AngleGrid.Size;
  for (auto &estimator : _Estimator)
    estimator.Initialize(
        {Para.Order + 1, Para.AngleGrid.Size, Para.KGrid.Size});
};

void ver4Obs::Measure0(double Factor) { Normalization += 1.0 * Factor; }

void ver4Obs::Measure(int Order, int QIndex, int AngleIndex,
                      const std::vector<verWeight> &Weight, double Factor) {

  ASSERT(Order != 0, "Order must be >=1!");

  ASSERT(AngleIndex >= 0 && AngleIndex < Para.AngleGrid.Size,
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
    VerFile << "# KGrid: " << Para.KGrid.ToString() << endl;
    VerFile << "# AngleGrid: " << Para.AngleGrid.ToString() << endl;

    for (int order = 0; order <= Para.Order; order++)
      for (int chan = 0; chan < 4; chan++)
        for (int angle = 0; angle < Para.AngleGrid.Size; ++angle)
          for (int qindex = 0; qindex < Para.KGrid.Size; ++qindex)
            for (int dir = 0; dir < 2; ++dir)
              VerFile << _Estimator[chan](order, angle, qindex)[dir] * PhyWeight
                      << "  ";
    VerFile.close();
  } else
    LOG_WARNING("Vertex4 for PID " << Para.PID << " fails to save!");
}
