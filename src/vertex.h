#ifndef vertex_H
#define vertex_H

// #include "utility/vector.h"
#include "global.h"
#include <Eigen/Dense>
#include <array>

double sum2(const momentum &);
double norm2(const momentum &);

namespace ver {
using Eigen::MatrixXd;

enum channel { I = 0, T, U, S };

const int MAXSIGMABIN = 100000;
const int SIGMA_TAU_BIN = 128;
const int SIGMA_K_BIN = 128;

class weightMatrix {
  // 2x2 matrix of weight; Direct/Exchange and chain/Lver/Rver/other
public:
  weightMatrix() { SetZero(); }
  void SetZero() {
    for (auto &i : _Weight)
      i = 0.0;
  }
  double Sum() {
    double sum = 0;
    for (auto &i : _Weight)
      sum += i;
    return sum;
  }

  double Abs() {
    double sum = 0;
    // for (auto &i : _Weight)
    //   sum += fabs(i);
    // sum = fabs(_Weight[DIR] + _Weight[EX] / 2.0) + fabs(_Weight[EX] / 2.0);
    sum = fabs(_Weight[DIR] + _Weight[EX] / SPIN);
    // sum = fabs(_Weight[DIR] + _Weight[EX]);
    return sum;
  }
  double &operator[](int dir) { return _Weight[dir]; }
  const double &operator[](int dir) const { return _Weight[dir]; }
  weightMatrix &operator+=(const weightMatrix &a) {
    _Weight[DIR] += a[DIR];
    _Weight[EX] += a[EX];
    return *this;
  }
  weightMatrix &operator*=(double Factor) {
    _Weight[DIR] *= Factor;
    _Weight[EX] *= Factor;
    return *this;
  }

  // weightMatrix operator*(double Factor) {
  //   weightMatrix Weight;
  //   Weight[DIR] = _Weight[DIR] *= Factor;
  //   Weight[EX] = _Weight[EX] *= Factor;
  //   return Weight;
  // }

private:
  array<double, 2> _Weight;
};

class fermi {
public:
  fermi();
  double Green(double Tau, const momentum &Momentum, spin Spin, int GType,
               double Scale = 0);

private:
  // beyond which the expense sigma function will be called
  double UpperBound, LowerBound;
  double DeltaK;
  double UpperBound2, LowerBound2; // lower upbound for better sigma
  double DeltaK2;
  double PhyGreen(double Tau, const momentum &Mom, int GType, double Scale = 0);
  double FockSigma(const momentum &Mom);
  double BuildFockSigma();
  double Fock(double k);
  // warning: this function only works for T=0 and 3D!!!!
  double GetSigma(double k);
  // double Sigma[MAXSIGMABIN];
  // double Sigma2[MAXSIGMABIN];

  // MatrixXd Sigma(SigmaMomBinSize, SigmaTauBinSize);
};

class verTensor {
public:
  verTensor();
  ~verTensor();
  void Initialize();
  double &Interaction(int Angle, int ExtQ, int Dir);
  double &Estimator(int Order, int Angle, int ExtQ, int Dir);

private:
  double *_Estimator;
  double *_Interaction;
  int QIndex;
  int AngleIndex;
  int OrderIndex;
};

// class sigmaTensor {
// public:
//   sigmaTensor();
//   ~sigmaTensor();
//   void Initialize();
//   double &Interaction(int Angle, int ExtQ, int Dir);
//   double &Estimator(int Order, int Angle, int ExtQ, int Dir);

// private:
// }

class verQTheta {
public:
  verQTheta();
  weightMatrix Interaction(const array<momentum *, 4> &LegK, double Tau,
                           bool IsRenorm, bool Boxed);

  void Measure(const momentum &InL, const momentum &InR, const int QIndex,
               int Order, const array<ver::weightMatrix, 4> &Weight,
               double Factor);
  void Update(double Ratio, int Order);
  void Save(bool Simple = false);
  void ClearStatis();
  void LoadWeight();
  void ResetIRScale(int IRScaleBin);

  array<verTensor, 4> Chan;

  // double TauBasis[TauBinSize][TauBasisNum];

  double Normalization;
  double PhyWeightT;
  double PhyWeightI;
};

double Angle3D(const momentum &K1, const momentum &K2);
double Index2Angle(const int &Index, const int &AngleNum);
int Angle2Index(const double &Angle, const int &AngleNum);
void _TestAngleIndex();
void _TestAngle2D();

double Index2Mom(const int &Index);
int Mom2Index(const double &K);

double Index2Scale(const int &Index);
int Scale2Index(const double &Scale);

double Index2Tau(const int &Index);
int Tau2Index(const double &Tau);

}; // namespace ver
#endif
