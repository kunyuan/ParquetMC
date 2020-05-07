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

  if (DiagType == POLAR) {
    PhyWeight = ExtMomBinSize * TauBinSize;
    Name = "polar";
  } else if (DiagType == SIGMA) {
    PhyWeight = ExtMomBinSize * TauBinSize;
    // PhyWeight = ExtMomBinSize / Para.Beta;
    Name = "sigma";
  } else if (DiagType == DELTA) {
    PhyWeight = ExtMomBinSize * TauBinSize;
    // PhyWeight = ExtMomBinSize;
    Name = "delta";
  } else
    return;

  _Estimator.Initialize({Para.Order + 1, ExtMomBinSize, TauBinSize});
}

void oneBodyObs::Measure0(double Factor) { Normalization += 1.0 * Factor; }
void oneBodyObs::Measure(int Order, int KBin, int TauBin, double Weight,
                         double Factor) {
  ASSERT(KBin >= 0 && KBin < ExtMomBinSize, "Kidx is out of range!");
  ASSERT(TauBin >= 0 && TauBin < TauBinSize, "TauIdx is out of range!");

  _Estimator(Order, KBin, TauBin) += Weight * Factor;
  _Estimator(0, KBin, TauBin) += Weight * Factor;
}

void oneBodyObs::Save() {

  for (int order = 0; order <= Para.Order; order++) {
    string FileName = fmt::format("{0}{1}_pid{2}.dat", Name, order, Para.PID);
    ofstream VerFile;
    VerFile.open(FileName, ios::out | ios::trunc);

    if (VerFile.is_open()) {

      VerFile << fmt::sprintf("#PID:%d, rs:%.3f, Beta: %.3f, Step: %d\n",
                              Para.PID, Para.Rs, Para.Beta, Var.Counter);

      VerFile << "# Norm: " << Normalization << endl;

      VerFile << "# TauTable: ";
      for (int t = 0; t < TauBinSize; ++t)
        VerFile << Para.Beta * t / TauBinSize << " ";

      VerFile << endl;
      VerFile << "# ExtMomBinTable: ";
      for (int qindex = 0; qindex < ExtMomBinSize; ++qindex)
        VerFile << Para.ExtMomTable[qindex][0] << " ";
      VerFile << endl;

      for (int tindex = 0; tindex < TauBinSize; ++tindex)
        for (int qindex = 0; qindex < ExtMomBinSize; ++qindex)
          VerFile << _Estimator(order, qindex, tindex) * PhyWeight << "  ";
      VerFile.close();
    } else {
      LOG_WARNING(Name << " for PID " << Para.PID << " fails to save!");
    }
  }
}

ver4Obs::ver4Obs() {
  Normalization = 1.0e-10;
  PhyWeight = AngBinSize;
  for (auto &estimator : _Estimator)
    estimator.Initialize({Para.Order + 1, AngBinSize, ExtMomBinSize});
};

void ver4Obs::Measure0(double Factor) { Normalization += 1.0 * Factor; }

void ver4Obs::Measure(const momentum &InL, const momentum &InR,
                      const int QIndex, int Order,
                      const std::vector<verWeight> &Weight, double Factor) {

  ASSERT(Order != 0, "Order must be >=1!");

  double CosAng = diag::Angle3D(InL, InR);
  int AngleIndex = diag::Angle2Index(CosAng, AngBinSize);

  ASSERT(AngleIndex >= 0 && AngleIndex < AngBinSize,
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
  for (int chan = 0; chan < 4; chan++) {
    for (int order = 0; order <= Para.Order; order++) {
      string FileName =
          fmt::format("vertex{0}_{1}_pid{2}.dat", order, chan, Para.PID);
      ofstream VerFile;
      VerFile.open(FileName, ios::out | ios::trunc);

      if (VerFile.is_open()) {

        VerFile << fmt::sprintf("#PID:%d, rs:%.3f, Beta: %.3f, Step: %d\n",
                                Para.PID, Para.Rs, Para.Beta, Var.Counter);

        VerFile << "# Norm: " << Normalization << endl;

        VerFile << "# AngleTable: ";
        for (int angle = 0; angle < AngBinSize; ++angle)
          VerFile << Para.AngleTable[angle] << " ";

        VerFile << endl;
        VerFile << "# ExtMomBinTable: ";
        for (int qindex = 0; qindex < ExtMomBinSize; ++qindex)
          VerFile << Para.ExtMomTable[qindex][0] << " ";
        VerFile << endl;

        for (int angle = 0; angle < AngBinSize; ++angle)
          for (int qindex = 0; qindex < ExtMomBinSize; ++qindex)
            for (int dir = 0; dir < 2; ++dir)
              VerFile << _Estimator[chan](order, angle, qindex)[dir] * PhyWeight
                      << "  ";
        VerFile.close();
      } else {
        LOG_WARNING("Vertex4 for PID " << Para.PID << " fails to save!");
      }
    }
  }
};