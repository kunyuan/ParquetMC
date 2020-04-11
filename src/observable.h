#ifndef observable_H
#define observable_H

#include "global.h"

namespace ver {
class sigData {
public:
  sigData(array<double, MaxTauNum> &TauTable);
  ~sigData();
  void Initialization();
  void Measure(int Order, const momentum &K, const vector<int> Tidx,
               const vector<double> Weight, double Factor); // all tau variables
  void Save();

private:
  const int TauNum = 128;
  const int KNum = 128;
  array<double, MaxTauNum> &Tau;
  double Normalization;
  double *_Estimator;
  double PhyWeight;

  int KIndex;
  int OrderIndex;
};

}; // namespace ver

#endif