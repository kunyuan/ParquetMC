#define FMT_HEADER_ONLY
#include "propagator.h"
#include "lib/green.h"
#include "utility/fmt/format.h"
#include "utility/fmt/printf.h"
#include <iostream>

using namespace diag;
using namespace std;
using namespace Eigen;

extern parameter Para;
extern variable Var;

void propagator::Initialize() {
  LoadInteraction();
  if (GreenType == BoldG) {
  }
}

double propagator::Green(double Tau, const momentum &K, spin Spin, int GType) {
  if (GType == 1)
    Tau = -1.0e-12;
  auto k = K.norm();
  auto Ek = k * k - Para.Mu; // bare propagator
  double fgreen = 0.0;

  if (GreenType == BoldG) {
    int BoldGOrder = Para.Order - 1;
    if (Para.Order == 2) {
      Ek += fockYukawa(k, Para.Kf, sqrt(Para.Lambda + Para.Mass2), true);
    } else if (Para.Order >= 3) {
      Ek +=
          _Interp1D<grid::FermiK>(_StaticSigma[BoldGOrder], _MomGridInterp, k);
    }
    fgreen = fermiGreen(Para.Beta, Tau, Ek);
    if (Para.Order >= 3)
      fgreen += _Interp2D<grid::FermiK>(_deltaGOrder[BoldGOrder],
                                        _MomGridInterp, k, Tau);
  } else if (GreenType == BareG) {
    fgreen = fermiGreen(Para.Beta, Tau, Ek);
  } else {
    Ek += fockYukawa(k, Para.Kf, sqrt(Para.Lambda + Para.Mass2), true);
    fgreen = fermiGreen(Para.Beta, Tau, Ek);
  }

  return fgreen;
}

void propagator::LoadGreen() {
  _DeltaG.setZero(Para.FermiKGrid.size, Para.TauGrid.size);

  // File.open("./selfconsistent/green.data", ios::in);
  // if (File.is_open()) {
  //   for (int k = 0; k < Para.FermiKGrid.size; ++k)
  //     for (int t = 0; t < Para.TauGrid.size; ++t)
  //       File >> _DeltaG(k, t);
  // } else {
  //   LOG_WARNING("Can not load Green weights! Initialze with zeros!\n");
  //   _DeltaG.setZero(Para.FermiKGrid.size, Para.TauGrid.size);
  // }
  // File.close();
  if (Para.Order > 2) {
    LoadGreenOrder();
    // SaveGreenOrder();
  }
}

void propagator::SaveGreenOrder() {
  ofstream File;
  std::array<momentum, 256> extKList;
  for (int k = 0; k < Para.FermiKGrid.size; ++k) {
    extKList[k].setZero();
    extKList[k][0] = Para.FermiKGrid.grid[k];
  }

  for (int o = 1; o <= Para.Order; o++) {
    char fname[100];
    snprintf(fname, 100, "SaveFullGreenOrder%d.data", o);
    File.open(fname, ios::out | ios::app);
    if (File.is_open()) {
      for (int k = 0; k < Para.FermiKGrid.size; ++k) {
        momentum K = extKList[k];
        for (int t = 0; t < Para.TauGrid.size - 1; ++t)
          File << Green(Para.TauGrid.grid[t], K, UP, 0) << "   ";
      }
    } else {
      LOG_WARNING("Fail to save the Green function's data.");
    }
    File.close();
  }
}

void propagator::LoadGreenOrder() {
  weight2D _deltaGTem;
  int _TauSize, _MomSize;
  double _MaxK;

  ifstream File1;
  File1.open("./selfconsistent/para.data", ios::in);
  if (File1.is_open()) {
    File1 >> _TauSize;
    File1 >> _MomSize;
    File1 >> _MaxK;
    _TauGridInterp.build(Para.Beta, _TauSize, 6.0 / Para.Ef);
    _MomGridInterp.build(Para.Kf, _MaxK * Para.Kf, _MomSize,
                         sqrt(1.0 / Para.Beta) * 2.0);
  } else {
    LOG_WARNING("Can not load para.data, FAILED\n");
    exit(0);
  }
  File1.close();

  weight1D _SigmaTem;
  ifstream File2;
  for (int o = 2; o <= Para.Order; o++) {
    _SigmaTem.setZero(_MomGridInterp.size);
    char fname[100];
    snprintf(fname, 100, "./selfconsistent/dispersion_order%d.data", o);
    File2.open(fname, ios::in);
    if (File2.is_open()) {
      for (int k = 0; k < _MomGridInterp.size; ++k)
        File2 >> _SigmaTem[k];
    } else {
      LOG_WARNING("Can not load dispersion! Initialze with zeros!\n");
      _SigmaTem.setZero(_MomGridInterp.size);
    }
    _StaticSigma[o] = _SigmaTem;
  }
  File2.close();

  ifstream File3;
  for (int o = 2; o <= Para.Order; o++) {
    _deltaGTem.setZero(_MomGridInterp.size, _TauGridInterp.size);
    char fname[100];
    snprintf(fname, 100, "./selfconsistent/green_order%d.data", o);
    File3.open(fname, ios::in);
    if (File3.is_open()) {
      for (int k = 0; k < _MomGridInterp.size; ++k) {
        for (int t = 0; t < _TauGridInterp.size; ++t)
          File3 >> _deltaGTem(k, t);
      }
    } else {
      LOG_WARNING("Can not load Order-Green weights! Initialze with zeros!\n");
      _deltaGTem.setZero(_MomGridInterp.size, _TauGridInterp.size);
    }
    _deltaGOrder[o] = _deltaGTem;
  }
  File3.close();
}

void propagator::LoadInteraction() {
  _DeltaRa.setZero(Para.BoseKGrid.size, Para.TauGrid.size);
  _DeltaRs.setZero(Para.BoseKGrid.size, Para.TauGrid.size);
  ifstream File;
  File.open("interaction.data", ios::in);
  if (File.is_open()) {
    for (int k = 0; k < Para.BoseKGrid.size; ++k)
      for (int t = 0; t < Para.TauGrid.size; ++t)
        File >> _DeltaRs(k, t);
    for (int k = 0; k < Para.BoseKGrid.size; ++k)
      for (int t = 0; t < Para.TauGrid.size; ++t)
        File >> _DeltaRa(k, t);
  } else {
    LOG_WARNING("Can not load interaction!");
  }
  // LOG_INFO("Getting: " << _DeltaRa(0, 0));
  // exit(0);
  File.close();
}

double propagator::F(double Tau, const momentum &K, spin Spin, int GType) {
  if (Tau == 0.0)
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

double lindhard(double x) {
  if (abs(x - 2) < 1.0e-5)
    return 0.5;
  else if (abs(x) < 1.0e-5)
    return 1.0;
  else
    return 0.5 - (x * x - 1) / 4.0 / x * log(abs((1 + x) / (1 - x)));
}

verWeight propagator::Interaction(const momentum &KInL, const momentum &KOutL,
                                  const momentum &KInR, const momentum &KOutR,
                                  double ExtQ) {
  verWeight Weight;
  // Weight = {1.0, 0.0};
  // return Weight;

  double kDiQ = (KInL - KOutL).norm();
  Weight[DIR] = -8.0 * PI * Para.Charge2 /
                (kDiQ * kDiQ + Para.Mass2 + Para.Lambda) / Para.Beta;

  if (DiagType == SIGMA && IsZero(kDiQ))
    Weight[DIR] = 0.0;

  // check irreducibility
  if ((IsProper || DiagType == POLAR) && IsEqual(kDiQ, ExtQ))
    Weight[DIR] = 0.0;

  double kExQ = (KInL - KOutR).norm();
  Weight[EX] = 8.0 * PI * Para.Charge2 /
               (kExQ * kExQ + Para.Mass2 + Para.Lambda) / Para.Beta;

  // Weight[EX] = 0.0;
  if (DiagType == SIGMA && IsZero(kExQ))
    Weight[EX] = 0.0;

  // check irreducibility
  if ((IsProper || DiagType == POLAR) && IsEqual(kExQ, ExtQ))
    Weight[EX] = 0.0;

  // cout << "Ver0: " << Weight[DIR] << ", " << Weight[EX] << endl;
  // cout << "extnal: " << ExtQ << ", " << kDiQ << endl;
  Weight[EX] = 0.0;
  return Weight;
}

verWeight propagator::InteractionTauBare(const momentum &KInL,
                                         const momentum &KOutL,
                                         const momentum &KInR,
                                         const momentum &KOutR, double inT,
                                         double outT, double ExtQ) {
  verWeight Weight;

  double kDiQ = (KInL - KOutL).norm();
  Weight[DIR] = -8.0 * PI /
                (kDiQ * kDiQ + Para.Mass2 +
                 Para.Nf * 8.0 * PI * lindhard(kDiQ / Para.Kf / 2.0)) /
                Para.Beta;

  Weight[DIR] += _Interp2D<grid::BoseK>(_DeltaRs, Para.BoseKGrid, Para.TauGrid,
                                        kDiQ, outT - inT);
  // Weight[DIR] =
  //     -8.0 * PI / (kDiQ * kDiQ + Para.Mass2 + Para.Nf * 8.0 * PI) /
  //     Para.Beta;

  double kExQ = (KInL - KOutR).norm();
  Weight[EX] = 8.0 * PI /
               (kExQ * kExQ + Para.Mass2 +
                Para.Nf * 8.0 * PI * lindhard(kDiQ / Para.Kf / 2.0));
  Weight[EX] -= _Interp2D<grid::BoseK>(_DeltaRs, Para.BoseKGrid, Para.TauGrid,
                                       kExQ, outT - inT);

  // Weight = {0.0, 0.0};
  // Weight[EX] = 0.0;
  return Weight;
}

verWeight propagator::InteractionTau(const momentum &KInL,
                                     const momentum &KOutL,
                                     const momentum &KInR,
                                     const momentum &KOutR, double inT,
                                     double outT, double ExtQ) {
  verWeight Weight;
  double Ws, Wa;

  double kDiQ = (KInL - KOutL).norm();
  Ws = _Interp2D<grid::BoseK>(_DeltaRs, Para.BoseKGrid, Para.TauGrid, kDiQ,
                              outT - inT);
  Wa = _Interp2D<grid::BoseK>(_DeltaRa, Para.BoseKGrid, Para.TauGrid, kDiQ,
                              outT - inT);
  Weight[DIR] = -(Ws - Wa);
  Weight[EX] = 2.0 * Wa;
  // cout << Ws << ", " << Wa << ", " << _DeltaRs(0, 0) << endl;

  double kExQ = (KInL - KOutR).norm();
  Ws = _Interp2D<grid::BoseK>(_DeltaRs, Para.BoseKGrid, Para.TauGrid, kExQ,
                              outT - inT);
  Wa = _Interp2D<grid::BoseK>(_DeltaRa, Para.BoseKGrid, Para.TauGrid, kExQ,
                              outT - inT);
  Weight[DIR] += -2.0 * Wa;
  Weight[EX] += Ws - Wa;

  // cout << "Ws: " << Ws << ", " << Wa << endl;
  // double v = 8.0 * PI / Para.Mass2;
  // Weight[DIR] = v * v * Para.Nf / (1.0 + v * Para.Nf) / Para.Beta;

  // Weight[EX] = 0.0;
  return Weight;
}

double propagator::Interaction(const momentum &TranQ, int VerOrder,
                               double ExtQ) {
  double kQ = TranQ.norm();
  if ((IsProper || DiagType == POLAR) && IsEqual(kQ, ExtQ))
    return 0.0;

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
      Weight *= pow(-1.0 * Weight * Para.Lambda / 8.0 / PI, VerOrder);
    return Weight;
  }
}

double propagator::CounterBubble(const momentum &K) {
  double Factor = Para.Lambda / (8.0 * PI * Para.Nf);
  // Factor *=
  //     Green(Para.Beta / 2.0, K, UP, 0) * Green(-Para.Beta / 2.0, K, UP, 0);

  double Ek = K.squaredNorm() - Para.Ef;
  Factor *= -0.5 / (1.0 + cosh(Ek * Para.Beta));
  // ASSERT_ALLWAYS(IsEqual())
  return Factor;
}

template <typename KGrid>
double propagator::_Interp1D(const weight1D &data, const KGrid &kgrid,
                             double K) {
  if (K > kgrid.grid.back() || K < kgrid.grid[0])
    return 0.0;

  int idx0 = kgrid.floor(K);
  int idx1 = idx0 + 1;
  double K0 = kgrid.grid[idx0];
  double K1 = kgrid.grid[idx1];

  return (data[idx0] * (K1 - K) + data[idx1] * (K - K0)) / (K1 - K0);
}

// template double progator::_Interp1D(const weight1D &data, double K,
//                                     grid::FermiK kgrid);

template <typename KGrid>
double propagator::_Interp2D(const weight2D &data, const KGrid &kgrid, double K,
                             double T) {
  double factor = 1.0;
  if (T < 0.0) {
    T = T + Para.Beta;
    factor *= -1.0;
  }

  if (K > kgrid.grid.back() || K < kgrid.grid[0] ||
      T > Para.TauGrid.grid.back() || T < 0.0)
    return 0.0;

  int Tidx0 = _TauGridInterp.floor(T);
  double dT0 = T - _TauGridInterp.grid[Tidx0],
         dT1 = _TauGridInterp.grid[Tidx0 + 1] - T;

  // ASSERT(dT0 >= 0.0 && dT1 >= 0.0,
  //  "Interpolate fails: " << T - dT0 << "<" << T << "<" << T + dT1);
  int Kidx0 = kgrid.floor(K);
  double dK0 = K - kgrid.grid[Kidx0], dK1 = kgrid.grid[Kidx0 + 1] - K;

  // ASSERT(dK0 >= 0.0 && dK1 >= 0.0,
  //  "Interpolate fails: " << K - dK0 << "<" << K << "<" << K + dK1);

  double d00 = data(Kidx0, Tidx0), d01 = data(Kidx0, Tidx0 + 1);
  double d10 = data(Kidx0 + 1, Tidx0), d11 = data(Kidx0 + 1, Tidx0 + 1);

  double g0 = d00 * dK1 + d10 * dK0;
  double g1 = d01 * dK1 + d11 * dK0;

  double gx = factor * (g0 * dT1 + g1 * dT0) / (dK0 + dK1) / (dT0 + dT1);
  return gx;
}

template <typename KGrid>
double propagator::_Interp2D(const weight2D &data, const KGrid &kgrid,
                             const grid::Tau &tgrid, double K, double T) {
  if (T < 0.0) {
    T = T + Para.Beta;
  }
  if (K > kgrid.grid.back() || T > Para.TauGrid.grid.back() || T < 0.0)
    return 0.0;

  int Tidx0 = tgrid.floor(T);
  double dT0 = T - tgrid.grid[Tidx0], dT1 = tgrid.grid[Tidx0 + 1] - T;

  if (K < kgrid.grid[0])
    return (data(0, Tidx0) * dT1 + data(0, Tidx0 + 1) * dT0) / (dT0 + dT1);

  // ASSERT(dT0 >= 0.0 && dT1 >= 0.0,
  //  "Interpolate fails: " << T - dT0 << "<" << T << "<" << T + dT1);
  int Kidx0 = kgrid.floor(K);
  double dK0 = K - kgrid.grid[Kidx0], dK1 = kgrid.grid[Kidx0 + 1] - K;

  // ASSERT(dK0 >= 0.0 && dK1 >= 0.0,
  //  "Interpolate fails: " << K - dK0 << "<" << K << "<" << K + dK1);

  double d00 = data(Kidx0, Tidx0), d01 = data(Kidx0, Tidx0 + 1);
  double d10 = data(Kidx0 + 1, Tidx0), d11 = data(Kidx0 + 1, Tidx0 + 1);

  double g0 = d00 * dK1 + d10 * dK0;
  double g1 = d01 * dK1 + d11 * dK0;

  // cout << K << "-" << Kidx0 << " dK0: " << dK0 << " dK1: " << dK1 << endl;
  // cout << d00 << ", " << d01 << ", " << d10 << ", " << d11 << endl;
  double gx = (g0 * dT1 + g1 * dT0) / (dK0 + dK1) / (dT0 + dT1);
  return gx;
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
