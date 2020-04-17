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
  _Estimator = new double[MaxOrder * TauNum * KNum];
  OrderIndex = TauNum * KNum;
  KIndex = TauNum;
  PhyWeight = 1.0;
  Initialization();
  return;
}

void sigData::Initialization() {
  for (int i = 0; i < MaxOrder * TauNum * KNum; ++i)
    _Estimator[i] = 0.0;
}

sigData::~sigData() { delete[] _Estimator; }

void sigData::Measure(int Order, const momentum &K, const vector<int> Tidx,
                      const vector<double> Weight, double Factor) {

  if (Order == 0) {
    Normalization += 1.0 * Factor;
    return;
  }
  int Size = Weight.size();
  ASSERT_ALLWAYS(Order == 0 && Size == 1, "Order 0 only has Size=1");
  int KIdx = int(K.norm() / Para.MaxExtMom * KNum);

  for (int i = 0; i < Size; ++i) {
    int TauIdx = int((Tau[Tidx[i]] / Para.Beta) * TauNum);
    _Estimator[Order * OrderIndex + KIdx * KIndex + TauIdx] +=
        Weight[i] * Factor;
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

      for (int qindex = 0; qindex < KNum; ++qindex)
        for (int tindex = 0; tindex < TauNum; ++tindex)
          VerFile << _Estimator[order * OrderIndex + qindex * KIndex + tindex] *
                         PhyWeight
                  << "  ";
      VerFile.close();
    } else {
      LOG_WARNING("Polarization for PID " << Para.PID << " fails to save!");
    }
  }
}