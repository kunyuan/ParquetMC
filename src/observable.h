#ifndef observable_H
#define observable_H

#include "global.h"

namespace ver {
class sigData {
public:
  sigData();
  ~sigData();
  void Initialization();
  void Measure0(double Factor);                          // all tau variables
  void Measure1(int Kidx, double Weight, double Factor); // all tau variables
  void Measure(int Order, int Kidx, const vector<int> Tidx,
               const vector<double> Weight, double Factor); // all tau variables
  void Save();
  void LoadWeight();

private:
  double Normalization;
  double *_Estimator;
  double *_EstimatorEqT;
  double *_EstimatorW1;
  double *_EstimatorW2;
  double PhyWeight;

  int KIndex;
  int OrderIndex;
};

class polarData {
public:
  polarData();
  ~polarData();
  void Initialization();
  void Measure(int Order, int Kidx, double Tau, double Weight,
               double Factor); // all tau variables
  void Save();

private:
  double Normalization;
  double *_Estimator;
  double PhyWeight;

  int KIndex;
  int OrderIndex;
};

class deltaData {
public:
  deltaData();
  ~deltaData();
  void Initialization();
  void Measure(int Order, int Kidx, double Tau, double Weight,
               double Factor); // all tau variables
  void Save();

private:
  double Normalization;
  double *_Estimator;
  double PhyWeight;

  int KIndex;
  int OrderIndex;
};

}; // namespace ver

#endif