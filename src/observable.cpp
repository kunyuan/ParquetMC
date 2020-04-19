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
  OrderIndex = TauBinSize * ExtMomBinSize;
  KIndex = TauBinSize;
  PhyWeight = 1.0;
  Initialization();
  return;
}

void sigData::Initialization() {
  for (int i = 0; i < MaxOrder * TauBinSize * ExtMomBinSize; ++i)
    _Estimator[i] = 0.0;
}

sigData::~sigData() { delete[] _Estimator; }

void sigData::Measure0(double Factor) { Normalization += 1.0 * Factor; }
void sigData::Measure1(int Kidx, double Weight, double Factor) {
  _Estimator[1 * OrderIndex + Kidx * KIndex + 0] += Weight * Factor;
  _Estimator[0 * OrderIndex + Kidx * KIndex + 0] += Weight * Factor;
}

void sigData::Measure(int Order, int Kidx, const vector<int> Tidx,
                      const vector<double> Weight, double Factor) {

  // only for Order >=2
  int Size = Weight.size();

  for (int i = 0; i < Size; ++i) {
    int TauIdx = int((Tau[Tidx[i]] / Para.Beta) * TauBinSize);
    _Estimator[Order * OrderIndex + Kidx * KIndex + TauIdx] +=
        Weight[i] * Factor;
    _Estimator[0 * OrderIndex + Kidx * KIndex + TauIdx] += Weight[i] * Factor;
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

      for (int qindex = 0; qindex < ExtMomBinSize; ++qindex)
        for (int tindex = 0; tindex < TauBinSize; ++tindex)
          VerFile << _Estimator[order * OrderIndex + qindex * KIndex + tindex] *
                         PhyWeight
                  << "  ";
      VerFile.close();
    } else {
      LOG_WARNING("Sigma for PID " << Para.PID << " fails to save!");
    }
  }
}

void sigData::LoadWeight() {}