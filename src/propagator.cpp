#define FMT_HEADER_ONLY
#include "propagator.h"
#include "utility/fmt/format.h"
#include "utility/fmt/printf.h"
#include <iostream>

using namespace diag;
using namespace std;
using namespace Eigen;

extern parameter Para;
extern variable Var;

double propagator::Green(double Tau, const momentum &K, spin Spin, int GType) {
  return _BareGreen(Tau, K, Spin, GType);
}

double propagator::_BareGreen(double Tau, const momentum &K, spin Spin,
                              int GType) {
  // if tau is exactly zero, set tau=0^-
  double green, Ek, kk, k;
  if (Tau == 0.0)
    Tau = -1.0e-12;

  // equal time green's function
  if (GType == 1)
    Tau = -1.0e-12;

  double s = 1.0;
  if (Tau < 0.0) {
    Tau += Para.Beta;
    s = -s;
  } else if (Tau >= Para.Beta) {
    Tau -= Para.Beta;
    s = -s;
  }

  k = K.norm();
  Ek = k * k; // bare propagator

  double x = Para.Beta * (Ek - Para.Mu) / 2.0;
  double y = 2.0 * Tau / Para.Beta - 1.0;
  if (x > 100.0)
    green = exp(-x * (y + 1.0));
  else if (x < -100.0)
    green = exp(x * (1.0 - y));
  else
    green = exp(-x * y) / (2.0 * cosh(x));

  green *= s;

  ASSERT(std::isnan(green) == false,
         "Step:" << Var.Counter << ", Green is too large! Tau=" << Tau
                 << ", Ek=" << Ek << ", Green=" << green << ", Mom=" << K);

  return green;
}

void propagator::LoadGreen(string FileName) {

  _StaticSigma.setZero(Para.KGrid.Size);
  _DeltaG.setZero(Para.KGrid.Size, Para.TauGrid.Size);

  ifstream File;
  File.open(FileName, ios::in);
  ASSERT_ALLWAYS(File.is_open(), "Can not load Green weights! \n");

  for (int k = 0; k < Para.KGrid.Size; ++k)
    File >> _StaticSigma[k];

  for (int k = 0; k < Para.KGrid.Size; ++k)
    for (int t = 0; t < Para.KGrid.Size; ++t)
      File >> _DeltaG(k, t);

  for (int k = 0; k < Para.KGrid.Size; ++k)
    cout << _StaticSigma[k] << ", " << endl;
}

double propagator::F(double Tau, const momentum &K, spin Spin, int GType) {
  if (Tau = 0.0)
    Tau = -1.0e-10;

  double Sign = -1.0;
  if (Tau < 0.0) {
    // make sure 0<Tau<Beta
    Tau = Para.Beta + Tau;
    Sign *= -1.0;
  }

  double Ek = K.squaredNorm() - Para.Mu;
  // return Sign * Tau * exp(-Ek * Tau) / 2.0 / (1 + cosh(Para.Beta * Ek));

  double Prefactor, Term1, Term2;
  if (Ek < 0.0) {
    double b = exp(Ek * Para.Beta);
    double x = exp(Ek * (Para.Beta - Tau));
    double y = 1.0 / (1 + 2.0 * b + b * b);
    Term1 = Tau * x * y;
    Term2 = 0.5 / Ek * exp(Ek * Tau) * (x * x - 1.0) * y;
    return (Term1 + Term2) * Sign;
  } else if (Ek > 0.0) {
    double b = exp(-Ek * Para.Beta);
    double x = exp(-Ek * (Para.Beta - Tau));
    double y = 1.0 / (1.0 + 2 * b + b * b);
    Term1 = Tau * exp(-Ek * (Tau + Para.Beta)) * y;
    Term2 = 0.5 / Ek * exp(-Ek * Tau) * (1.0 - x * x) * y;
    return (Term1 + Term2) * Sign;
  } else
    // unphysical
    return 0.0;
}

verWeight propagator::Interaction(const momentum &KInL, const momentum &KOutL,
                                  const momentum &KInR, const momentum &KOutR,
                                  bool Boxed, double ExtQ) {
  verWeight Weight;
  // Weight = {1.0, 0.0};
  // return Weight;

  double kDiQ = (KInL - KOutL).norm();
  Weight[DIR] =
      -8.0 * PI * Para.Charge2 / (kDiQ * kDiQ + Para.Mass2 + Para.Lambda);

  if (DiagType == SIGMA && IsZero(kDiQ))
    Weight[DIR] = 0.0;

  // check irreducibility
  if (DiagType == POLAR && IsEqual(kDiQ, ExtQ))
    Weight[DIR] = 0.0;

  if (!Boxed) {
    double kExQ = (KInL - KOutR).norm();
    Weight[EX] =
        8.0 * PI * Para.Charge2 / (kExQ * kExQ + Para.Mass2 + Para.Lambda);

    if (DiagType == SIGMA && IsZero(kExQ))
      Weight[EX] = 0.0;

    // check irreducibility
    if (DiagType == POLAR && IsEqual(kExQ, ExtQ))
      Weight[EX] = 0.0;

  } else
    Weight[EX] = 0.0;

  // cout << "Ver0: " << Weight[DIR] << ", " << Weight[EX] << endl;
  // cout << "extnal: " << ExtQ << ", " << kDiQ << endl;
  return Weight;
}

double propagator::Interaction(const momentum &TranQ, int VerOrder) {
  double kQ = TranQ.norm();
  if (VerOrder < 0) {
    // Bare interaction
    if (kQ > 1.0e-8)
      return -8.0 * PI * Para.Charge2 / (kQ * kQ);
    else
      return 0.0;
  } else {
    // Order N shifted interaction
    double Weight =
        -8.0 * PI * Para.Charge2 / (kQ * kQ + Para.Mass2 + Para.Lambda);
    if (VerOrder > 0)
      Weight *= pow(Weight * Para.Lambda / 8.0 / PI, VerOrder);
    return Weight;
  }
}

double propagator::CounterBubble(const momentum &K) {
  double Factor = Para.Lambda / (8.0 * PI * Para.Nf);
  Factor *=
      Green(Para.Beta / 2.0, K, UP, 0) * Green(-Para.Beta / 2.0, K, UP, 0);
  return Factor;
}

double propagator::_Interp1D(double K, const weight1D &data) {
  int idx0 = Para.KGrid.Floor(K);
  int idx1 = idx0 + 1;
  double K0 = Para.KGrid.Grid[idx0];
  double K1 = Para.KGrid.Grid[idx1];
  ASSERT(K0 <= K && K1 > K,
         "Interpolate fails: " << K0 << "<" << K << "<" << K1);
  return (data[idx0] * (K1 - K) + data[idx1] * (K - K0)) / (K1 - K0);
}
double propagator::_Interp2D(double K, double T, const weight2D &data) {
  int Tidx0 = Para.TauGrid.Floor(T);
  double dT0 = T - Para.TauGrid.Grid[Tidx0],
         dT1 = Para.TauGrid.Grid[Tidx0 + 1] - T;
  ASSERT(dT0 >= 0.0 && dT1 > 0.0,
         "Interpolate fails: " << T - dT0 << "<" << T << "<" << T + dT1);

  int Kidx0 = Para.KGrid.Floor(K);
  int dK0 = K - Para.KGrid.Grid[Kidx0], dK1 = Para.KGrid.Grid[Kidx0 + 1] - K;
  ASSERT(dK0 >= 0.0 && dK1 > 0.0,
         "Interpolate fails: " << K - dK0 << "<" << K << "<" << K + dK1);

  double d00 = data(Kidx0, Tidx0), d01 = data(Kidx0, Tidx0 + 1);
  double d10 = data(Kidx0 + 1, Tidx0), d11 = data(Kidx0 + 1, Tidx0 + 1);

  double g0 = d00 * dK1 + d10 * dK0;
  double g1 = d01 * dK1 + d11 * dK0;
  return (g0 * dT1 + g1 * dT0) / (dK0 + dK1) / (dT0 + dT1);
}

double diag::Angle3D(const momentum &K1, const momentum &K2) {
  // Returns the angle in radians between vectors 'K1' and 'K2'
  double dotp = K1.dot(K2);
  double Angle2D = dotp / K1.norm() / K2.norm();
  return Angle2D;
}

void diag::_TestAngle2D() {
  // Test Angle functions
  momentum K1, K2;
  K1.setZero();
  K2.setZero();
  K1[0] = 1.0;
  K2[0] = 1.0;

  ASSERT_ALLWAYS(
      abs(Angle3D(K1, K2) - 1.0) < 1.e-7,
      fmt::format("Angle between K1 and K2 are not zero! It is {:.13f}",
                  Angle3D(K1, K2)));

  K1.setZero();
  K2.setZero();
  K1[0] = 1.0;
  K2[0] = -1.0;
  ASSERT_ALLWAYS(
      abs(Angle3D(K1, K2) - (-1.0)) < 1.e-7,
      fmt::format("Angle between K1 and K2 are not Pi! Instead, it is {:.13f}",
                  Angle3D(K1, K2)));

  K1.setZero();
  K2.setZero();
  K1[0] = 1.0;
  K2[0] = -1.0;
  K1[1] = 0.0;
  K2[1] = -EPS;
  ASSERT_ALLWAYS(
      abs(Angle3D(K1, K2) - 1.0) < 1.e-7,
      fmt::format("Angle between K1 and K2 are not 2.0*Pi! It is {:.13f}",
                  Angle3D(K1, K2)));
}