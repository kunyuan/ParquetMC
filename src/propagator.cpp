#include "propagator.h"

using namespace diag;
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
    Tau = -1.0e-10;

  // equal time green's function
  if (GType == 1)
    Tau = -1.0e-10;

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

verWeight propagator::Interaction(const momentum &KInL, const momentum &KOutL,
                                  const momentum &KInR, const momentum &KOutR,
                                  bool Boxed, double ExtQ) {

  verWeight Weight;

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
  double Weight =
      -8.0 * PI * Para.Charge2 / (kQ * kQ + Para.Mass2 + Para.Lambda);

  // interactions with counterterms
  if (VerOrder > 0)
    Weight *= pow(Weight * Para.Lambda / 8.0 / PI, VerOrder);

  return Weight;
}