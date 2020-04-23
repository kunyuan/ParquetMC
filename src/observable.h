#ifndef observable_H
#define observable_H

#include "global.h"

namespace ver {
class sigData {
public:
  sigData(array<double, MaxTauNum> &TauTable);
  ~sigData();
  void Initialization();
  void Measure0(double Factor);                          // all tau variables
  void Measure1(int Kidx, double Weight, double Factor); // all tau variables
  void Measure(int Order, int Kidx, const vector<int> Tidx,
               const vector<double> Weight, double Factor); // all tau variables
  void Save();
  void LoadWeight();

private:
  array<double, MaxTauNum> &Tau;
  double Normalization;
  double *_Estimator;
  double *_EstimatorEqT;
  double *_EstimatorW1;
  double *_EstimatorW2;
  double PhyWeight;

  int KIndex;
  int OrderIndex;
};

}; // namespace ver

#endif