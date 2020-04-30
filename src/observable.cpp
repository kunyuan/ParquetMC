#include "observable.h"
#include "global.h"
#include "utility/abort.h"
#include "utility/fmt/format.h"
#include "utility/fmt/printf.h"
#include "utility/utility.h"
#include <cmath>
#include <iostream>

using namespace ver;

extern parameter Para;

sigData::sigData(array<double, MaxTauNum> &TauTable) : Tau(TauTable) {
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
      int TauIdx = int((Tau[T[i]] / Para.Beta) * TauBinSize);
      _Estimator[Order * OrderIndex + Kidx * KIndex + TauIdx] +=
          Weight[i] * Factor * TauBinSize / Para.Beta;
      _Estimator[0 * OrderIndex + Kidx * KIndex + TauIdx] +=
          Weight[i] * Factor * TauBinSize / Para.Beta;

      _EstimatorW1[Order * ExtMomBinSize + Kidx] +=
          Weight[i] * Factor * sin(PI * (Tau[T[i]] - Tau[0]) / Para.Beta);
      _EstimatorW2[Order * ExtMomBinSize + Kidx] +=
          Weight[i] * Factor * sin(-PI * (Tau[T[i]] - Tau[0]) / Para.Beta);
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
                              Para.PID, Para.Rs, Para.Beta, Para.Counter);

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

polarData::polarData() {
  _Estimator = new double[MaxOrder * TauBinSize * ExtMomBinSize];
  OrderIndex = TauBinSize * ExtMomBinSize;
  KIndex = TauBinSize;
  PhyWeight = ExtMomBinSize;
  Initialization();
  return;
}

void polarData::Initialization() {
  for (int i = 0; i < MaxOrder * TauBinSize * ExtMomBinSize; ++i)
    _Estimator[i] = 0.0;
}

polarData::~polarData() { delete[] _Estimator; }

void polarData::Measure(int Order, int Kidx, double Tau, double Weight,
                        double Factor) {

  if (Order == 0)
    Normalization += 1.0 * Factor;
  // only for Order >=1
  else {
    if (Tau < 0.0)
      Tau = Tau + Para.Beta;

    int TauIdx = int(Tau / Para.Beta * TauBinSize);

    _Estimator[Order * OrderIndex + Kidx * KIndex + TauIdx] +=
        Weight * Factor * TauBinSize / Para.Beta;
    _Estimator[0 * OrderIndex + Kidx * KIndex + TauIdx] +=
        Weight * Factor * TauBinSize / Para.Beta;
  }
}

void polarData::Save() {

  for (int order = 0; order <= Para.Order; order++) {
    string FileName = fmt::format("polar{0}_pid{1}.dat", order, Para.PID);
    ofstream VerFile;
    VerFile.open(FileName, ios::out | ios::trunc);

    if (VerFile.is_open()) {

      VerFile << fmt::sprintf("#PID:%d, rs:%.3f, Beta: %.3f, Step: %d\n",
                              Para.PID, Para.Rs, Para.Beta, Para.Counter);

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
      LOG_WARNING("Polar for PID " << Para.PID << " fails to save!");
    }
  }
}

void polarData::LoadWeight() {}