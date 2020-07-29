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
  void Measure(int Order, int QIndex, int AngIdx,
               const std::vector<verWeight> &Weight, double Factor);
  void Reset();
  void Save();

private:
  double Normalization;
  std::array<tensor3<verWeight>, 4> _Estimator;
  double PhyWeight;
};

class oneBodyObs {
public:
  oneBodyObs();
  void Measure0(double Factor);
  void Measure(int Order, int KBin, int TauBin, double Weight,
               double Factor); // all tau variables
  void Reset();
  void Save();
  void Save(int);
  void Check();
private:
  std::string Name;
  int channel; // channel s=0,p=1,d=2 ...
  double Normalization;
  tensor3<double> _Estimator;
  double PhyWeight;
  int ksize;
};

} // namespace obs

#endif
