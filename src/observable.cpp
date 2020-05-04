#include "observable.h"
#include "global.h"
#include "propagator.h"
#include "utility/abort.h"
#include "utility/fmt/format.h"
#include "utility/fmt/printf.h"
#include "utility/utility.h"
#include <cmath>
#include <iostream>

using namespace ver;
using namespace std;
using namespace obs;

extern parameter Para;
extern variable Var;

sigData::sigData() {
  _Estimator = new double[MaxOrder * TauBinSize * ExtMomBinSize];
  _EstimatorEqT = new double[MaxOrder * ExtMomBinSize];
  _EstimatorW1 = new double[MaxOrder * ExtMomBinSize];
  _EstimatorW2 = new double[MaxOrder * ExtMomBinSize];
  OrderIndex = TauBinSize * ExtMomBinSize;
  KIndex = TauBinSize;
  PhyWeight = ExtMomBinSize /
              Para.Beta; // order 1 sigma doesn't have a real tau variable
  Initialization();
  return;
}

void sigData::Initialization() {
  for (int i = 0; i < MaxOrder * TauBinSize * ExtMomBinSize; ++i)
    _Estimator[i] = 0.0;
  for (int i = 0; i < MaxOrder * ExtMomBinSize; ++i) {
    _EstimatorEqT[i] = 0.0;
    _EstimatorW1[i] = 0.0;
    _EstimatorW2[i] = 0.0;
  }
}

sigData::~sigData() {
  delete[] _Estimator;
  delete[] _EstimatorEqT;
  delete[] _EstimatorW1;
  delete[] _EstimatorW2;
}

void sigData::Measure0(double Factor) { Normalization += 1.0 * Factor; }
void sigData::Measure1(int Kidx, double Weight, double Factor) {
  // cout << Kidx << endl;
  _EstimatorEqT[1 * ExtMomBinSize + Kidx] += Weight * Factor;
  _EstimatorEqT[0 * ExtMomBinSize + Kidx] += Weight * Factor;
}

void sigData::Measure(int Order, int Kidx, const vector<int> T,
                      const vector<double> Weight, double Factor) {

  // only for Order >=2
  int Size = T.size();

  for (int i = 0; i < Size; ++i) {
    // cout << Tidx[i] << endl;
    // cout << "i=" << i << ", " << T[i] << ", " << Tau[T[i]] << ", " <<
    // Weight[i]
    //      << endl;
    if (T[i] == 0) {
      _EstimatorEqT[Order * ExtMomBinSize + Kidx] += Weight[i] * Factor;
      _EstimatorEqT[0 * ExtMomBinSize + Kidx] += Weight[i] * Factor;
    } else {
      int TauIdx =
          int(((Var.Tau[T[i]] - Var.Tau[T[0]]) / Para.Beta) * TauBinSize);
      _Estimator[Order * OrderIndex + Kidx * KIndex + TauIdx] +=
          Weight[i] * Factor * TauBinSize / Para.Beta;
      _Estimator[0 * OrderIndex + Kidx * KIndex + TauIdx] +=
          Weight[i] * Factor * TauBinSize / Para.Beta;

      _EstimatorW1[Order * ExtMomBinSize + Kidx] +=
          Weight[i] * Factor *
          sin(PI * (Var.Tau[T[i]] - Var.Tau[0]) / Para.Beta);
      _EstimatorW2[Order * ExtMomBinSize + Kidx] +=
          Weight[i] * Factor *
          sin(-PI * (Var.Tau[T[i]] - Var.Tau[0]) / Para.Beta);
    }
  }
}

void sigData::Save() {

  for (int order = 0; order <= Para.Order; order++) {
    string FileName = fmt::format("sigma{0}_pid{1}.dat", order, Para.PID);
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

      for (int qindex = 0; qindex < ExtMomBinSize; ++qindex)
        VerFile << _EstimatorEqT[order * ExtMomBinSize + qindex] * PhyWeight
                << "  ";

      for (int qindex = 0; qindex < ExtMomBinSize; ++qindex)
        VerFile << _EstimatorW1[order * ExtMomBinSize + qindex] * PhyWeight
                << "  ";

      for (int qindex = 0; qindex < ExtMomBinSize; ++qindex)
        VerFile << _EstimatorW2[order * ExtMomBinSize + qindex] * PhyWeight
                << "  ";

      for (int tindex = 0; tindex < TauBinSize; ++tindex)
        for (int qindex = 0; qindex < ExtMomBinSize; ++qindex)
          VerFile << _Estimator[order * OrderIndex + qindex * KIndex + tindex] *
                         PhyWeight
                  << "  ";
      // VerFile << 0.0 << "  ";
      VerFile.close();
    } else {
      LOG_WARNING("Sigma for PID " << Para.PID << " fails to save!");
    }
  }
}

void sigData::LoadWeight() {}

polarObs::polarObs() {
  Normalization = 1.0e-10;
  PhyWeight = TauBinSize;
  _Estimator.Initialize({Para.Order + 1, TauBinSize, ExtMomBinSize});
}

void polarObs::Measure0(double Factor) { Normalization += 1.0 * Factor; }
void polarObs::Measure(int Order, int Kidx, double Tau, double Weight,
                       double Factor) {

  if (Order == 0)
    Normalization += 1.0 * Factor;
  // only for Order >=1
  else {
    if (Tau < 0.0)
      Tau = Tau + Para.Beta;

    int TauIdx = int(Tau / Para.Beta * TauBinSize);

    ASSERT(Kidx >= 0 && Kidx < ExtMomBinSize, "Kidx is out of range!");
    ASSERT(TauIdx >= 0 && TauIdx < TauBinSize, "TauIdx is out of range!");

    _Estimator(Order, Kidx, TauIdx) += Weight * Factor * TauBinSize / Para.Beta;
    _Estimator(0, Kidx, TauIdx) += Weight * Factor * TauBinSize / Para.Beta;
  }
}

void polarObs::Save() {

  for (int order = 0; order <= Para.Order; order++) {
    string FileName = fmt::format("polar{0}_pid{1}.dat", order, Para.PID);
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
      LOG_WARNING("Polar for PID " << Para.PID << " fails to save!");
    }
  }
}

deltaData::deltaData() {
  _Estimator = new double[MaxOrder * TauBinSize * ExtMomBinSize];
  OrderIndex = TauBinSize * ExtMomBinSize;
  KIndex = TauBinSize;
  PhyWeight = ExtMomBinSize / Para.Beta;
  Initialization();
  return;
}

void deltaData::Initialization() {
  for (int i = 0; i < MaxOrder * TauBinSize * ExtMomBinSize; ++i)
    _Estimator[i] = 0.0;
}

deltaData::~deltaData() { delete[] _Estimator; }

void deltaData::Measure(int Order, int Kidx, double Tau, double Weight,
                        double Factor) {

  if (Order == 0)
    Normalization += 1.0 * Factor;
  // only for Order >=1
  else {
    int TauIdx = 0;
    // if Order==1, then F must a delta function of tau

    if (Order >= 2) {
      if (Tau < 0.0)
        Tau = Tau + Para.Beta;
      TauIdx = int(Tau / Para.Beta * TauBinSize);
    }

    _Estimator[Order * OrderIndex + Kidx * KIndex + TauIdx] +=
        Weight * Factor * TauBinSize / Para.Beta;
    _Estimator[0 * OrderIndex + Kidx * KIndex + TauIdx] +=
        Weight * Factor * TauBinSize / Para.Beta;
  }
}

void deltaData::Save() {

  for (int order = 0; order <= Para.Order; order++) {
    string FileName = fmt::format("delta{0}_pid{1}.dat", order, Para.PID);
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
          VerFile << _Estimator[order * OrderIndex + qindex * KIndex + tindex] *
                         PhyWeight
                  << "  ";
      // VerFile << 0.0 << "  ";
      VerFile.close();
    } else {
      LOG_WARNING("Delta for PID " << Para.PID << " fails to save!");
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