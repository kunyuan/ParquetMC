#include "global.h"
#include "utility/abort.h"
#include "utility/fmt/format.h"
#include "utility/timer.h"
#include "weight.h"
#include <array>
#include <iostream>
#include <string>
#include <vector>

using namespace diag;
using namespace std;

extern parameter Para;
extern variable Var;
extern propagator Prop;

void weight::Test(int LoopNum) {
  // cout << "start testing ..." << endl;
  if (DiagType == GAMMA && LoopNum == 1)
    _TestOneLoopGamma();
  else if (DiagType == SIGMA && LoopNum == 2)
    _TestTwoLoopSigma();
}

void weight::_TestOneLoopGamma() {
  momentum InL = Var.LoopMom[INL];
  momentum OutL = Var.LoopMom[OUTL];
  momentum InR = Var.LoopMom[INR];
  momentum OutR = Var.LoopMom[OUTR];
  momentum K0 = Var.LoopMom[4];
  double Factor = 1.0 / pow(2.0 * PI, D);
  double dTau, GWeight, GWeightInBox;
  verWeight TestWeight, Weight;

  // one loop T and TC diagram
  momentum K1 = OutL + K0 - InL;
  auto LVer = Prop.Interaction(InL, OutL, K1, K0, false);
  auto RVer = Prop.Interaction(K0, K1, InR, OutR, false);
  dTau = Var.Tau[1] - Var.Tau[0];
  GWeight = Prop.Green(dTau, K0, UP, 0) * Prop.Green(-dTau, K1, UP, 0);
  GWeightInBox = Prop.Green(dTau, K0, UP, 0) * Prop.Green(-dTau, K0, UP, 0);

  TestWeight.setZero();
  TestWeight[DIR] = LVer[DIR] * RVer[DIR] * SPIN + LVer[EX] * RVer[DIR] +
                    LVer[DIR] * RVer[EX];
  TestWeight[EX] = LVer[EX] * RVer[EX];
  TestWeight *= GWeight;
  TestWeight[DIR] -= GWeightInBox * LVer[DIR] * RVer[DIR] * Para.Lambda /
                     (8.0 * PI) / Para.Nf * SPIN;
  TestWeight *= Factor * SymFactor[T];

  Weight = _GetWeight(1, {T, TC});

  if (abs(Weight[DIR] - TestWeight[DIR]) > 1.0e-10 ||
      abs(Weight[EX] - TestWeight[EX]) > 1.0e-10) {
    cout << fmt::format("G: {}", GWeight * Factor) << endl;
    cout << fmt::format(
                "GInBox: {}, sep: {}, {}",
                GWeightInBox * Factor * Para.Lambda / (8.0 * PI) / Para.Nf,
                Prop.Green(dTau, K0, UP, 0), Prop.Green(-dTau, K0, UP, 0))
         << endl;
    cout << fmt::format("LVer: {}, {}", LVer[DIR], LVer[EX]) << endl;
    cout << fmt::format("RVer: {}, {}", RVer[DIR], RVer[EX]) << endl;
    ABORT(fmt::format(
        "T channel one loop Weight Calculation has a bug! ({},{}) vs ({}, {})",
        Weight[DIR], Weight[EX], TestWeight[DIR], TestWeight[EX]));
  }

  // one loop U and UC diagram
  momentum K2 = OutR + K0 - InL;
  LVer = Prop.Interaction(InL, OutR, K2, K0, false);
  RVer = Prop.Interaction(K0, K2, InR, OutL, false);
  dTau = Var.Tau[1] - Var.Tau[0];
  GWeight = Prop.Green(dTau, K0, UP, 0) * Prop.Green(-dTau, K2, UP, 0);
  GWeightInBox = Prop.Green(dTau, K0, UP, 0) * Prop.Green(-dTau, K0, UP, 0);

  TestWeight.setZero();
  TestWeight[DIR] = LVer[EX] * RVer[EX];
  TestWeight[EX] = LVer[DIR] * RVer[DIR] * SPIN + LVer[EX] * RVer[DIR] +
                   LVer[DIR] * RVer[EX];
  TestWeight *= GWeight;
  TestWeight[EX] -= GWeightInBox * LVer[DIR] * RVer[DIR] * Para.Lambda /
                    (8.0 * PI) / Para.Nf * SPIN;
  TestWeight *= Factor * SymFactor[U];

  Weight = _GetWeight(1, {U, UC});
  if (abs(Weight[DIR] - TestWeight[DIR]) > 1.0e-10 ||
      abs(Weight[EX] - TestWeight[EX]) > 1.0e-10) {
    cout << fmt::format("G: {}", GWeight * Factor) << endl;
    cout << fmt::format(
                "GInBox: {}, sep: {}, {}",
                GWeightInBox * Factor * Para.Lambda / (8.0 * PI) / Para.Nf,
                Prop.Green(dTau, K0, UP, 0), Prop.Green(-dTau, K0, UP, 0))
         << endl;
    cout << fmt::format("LVer: {}, {}", LVer[DIR], LVer[EX]) << endl;
    cout << fmt::format("RVer: {}, {}", RVer[DIR], RVer[EX]) << endl;
    ABORT(fmt::format(
        "U channel one loop Weight Calculation has a bug! ({},{}) vs ({}, {})",
        Weight[DIR], Weight[EX], TestWeight[DIR], TestWeight[EX]));
  }

  // one loop S diagram
  momentum K3 = InR + InL - K0;
  LVer = Prop.Interaction(InL, K3, InR, K0, false);
  RVer = Prop.Interaction(K0, OutL, K3, OutR, false);
  dTau = Var.Tau[1] - Var.Tau[0];
  GWeight = Prop.Green(dTau, K0, UP, 0) * Prop.Green(dTau, K3, UP, 0);

  TestWeight.setZero();
  TestWeight[DIR] = LVer[EX] * RVer[DIR] + LVer[DIR] * RVer[EX];
  TestWeight[EX] = LVer[DIR] * RVer[DIR] + LVer[EX] * RVer[EX];
  TestWeight *= GWeight;
  TestWeight *= Factor * SymFactor[S];

  //   cout << "Calculate S chanel" << endl;
  Weight = _GetWeight(1, {S});
  if (abs(Weight[DIR] - TestWeight[DIR]) > 1.0e-10 ||
      abs(Weight[EX] - TestWeight[EX]) > 1.0e-10) {
    cout << fmt::format("G: {} G0 {} G3 {}", GWeight * Factor * SymFactor[S],
                        Prop.Green(dTau, K0, UP, 0),
                        Prop.Green(dTau, K3, UP, 0))
         << endl;
    cout << fmt::format("LVer: {}, {}", LVer[DIR], LVer[EX]) << endl;
    cout << fmt::format("RVer: {}, {}", RVer[DIR], RVer[EX]) << endl;
    ABORT(fmt::format(
        "S channel one loop Weight Calculation has a bug! ({},{}) vs ({}, {})",
        Weight[DIR], Weight[EX], TestWeight[DIR], TestWeight[EX]));
  }

  // cout << "Test End" << endl;
  return;
}

void weight::_TestTwoLoopGamma() {
  // momentum InL = Var.LoopMom[INL];
  // momentum OutL = Var.LoopMom[OUTL];
  // momentum InR = Var.LoopMom[INR];
  // momentum OutR = Var.LoopMom[OUTR];
  // momentum K0 = Var.LoopMom[4];
  // double Factor = 1.0 / pow(2.0 * PI, D);
  // double dTau, GWeight, GWeightInBox;
  // ver::weightMatrix TestWeight, Weight;

  // return;
}

double weight::_TestTwoLoopSigma() {
  double Factor = 1.0 / pow(2.0 * PI, D);
  momentum &ExtK = Var.LoopMom[0];
  momentum &K1 = Var.LoopMom[1];
  momentum &K2 = Var.LoopMom[2];
  double G1 = Prop.Green(Var.Tau[1] - Var.Tau[0], K1, UP, 0);
  double G2 = Prop.Green(Var.Tau[1] - Var.Tau[0], K2, UP, 0);
  double G3 = Prop.Green(Var.Tau[0] - Var.Tau[1], K1 + K2 - ExtK, UP, 0);
  double VerWeightDir = Prop.Interaction(K1 - ExtK);
  double VerWeightExLeft = Prop.Interaction(K2 - ExtK);

  double Weight1 =
      -VerWeightDir * VerWeightDir * G1 * G2 * G3 * SPIN * Factor * Factor;
  double Weight2 =
      VerWeightDir * VerWeightExLeft * G1 * G2 * G3 * Factor * Factor;
  double Weight3 =
      VerWeightExLeft * VerWeightExLeft * G1 * G2 * G3 * SPIN * Factor * Factor;
  // cout << "G_test=" << G1 << ", " << G2 << ", " << G3 << endl;
  // cout << "Ver_test=" << VerWeightDir << ", " << VerWeightExLeft << endl;

  cout << "DIR=" << Weight1 + 2 * Weight2 << ", EX=" << -Weight3 << endl;
  return Weight1 + 2 * Weight2;
}

verWeight weight::_GetWeight(int LoopNum, std::vector<channel> Channel) {
  verWeight Weight;
  // array<momentum *, 4> ExtLegK = {&Var.LoopMom[0], &Var.LoopMom[1],
  //                                 &Var.LoopMom[2], &Var.LoopMom[3]};
  vertex4 Ver4;
  Ver4.Build(0, // level
             1, // loopNum
             4, // loop index of the first internal K
             0, // tau index of the InTL leg
             Channel, RIGHT, false);

  Ver4.Evaluate(Var.LoopMom[INL], Var.LoopMom[OUTL], Var.LoopMom[INR],
                Var.LoopMom[OUTR], true);
  auto &ChanWeight = Ver4.ChanWeight;
  for (auto &w : ChanWeight)
    // collapse all channel to I
    ChanWeight[0] += w;
  // cout << ChanWeight[0].Abs() << endl;
  Weight = ChanWeight[0];
  return Weight;
}
