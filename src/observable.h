#ifndef observable_H
#define observable_H

#include "global.h"
#include <vector>

namespace obs {

enum channel { I = 0, T, U, S };

template <typename T> class tensor3 {
public:
  tensor3() { _Estimator = nullptr; };
  ~tensor3() {
    if (_Estimator != nullptr)
      delete[] _Estimator;
  };
  void Initialize(std::array<int, 3> dim) {
    Size = dim[0] * dim[1] * dim[2];
    Index0 = dim[1] * dim[2];
    Index1 = dim[2];
    _Estimator = new T[Size];
    for (int i = 0; i < Size; ++i)
      _Estimator[i] *= 0.0; // set all elements to zero
  };
  T &operator()(int X, int Y, int Z) {
    ASSERT(X * Index0 + Y * Index1 + Z < Size,
           "Estimator out of range! " << X << "x" << Y << "x" << Z);
    return _Estimator[X * Index0 + Y * Index1 + Z];
  };

private:
  T *_Estimator;
  int Index0, Index1, Size;
};

class ver4Obs {
public:
  ver4Obs();
  void Measure0(double Factor);
  void Measure(const momentum &InL, const momentum &InR, const int QIndex,
               int Order, const std::vector<verWeight> &Weight, double Factor);
  void Save();

private:
  double Normalization;
  std::array<tensor3<verWeight>, 4> _Estimator;
  double PhyWeight;
};
} // namespace obs

namespace ver {
class sigData {
public:
  sigData();
  ~sigData();
  void Initialization();
  void Measure0(double Factor);                          // all tau variables
  void Measure1(int Kidx, double Weight, double Factor); // all tau variables
  void Measure(int Order, int Kidx, const std::vector<int> Tidx,
               const std::vector<double> Weight,
               double Factor); // all tau variables
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