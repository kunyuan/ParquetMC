#define FMT_HEADER_ONLY
#include "utility/fmt/format.h"
#include "vertex4.h"
#include <array>
#include <iostream>
#include <string>
#include <vector>

using namespace diag;
using namespace std;

extern parameter Para;
extern variable Var;
extern propagator Prop;

void vertex4::_TestOneLoopGamma() {
  momentum InL = Var.LoopMom[INL];
  momentum OutL = Var.LoopMom[OUTL];
  momentum InR = Var.LoopMom[INR];
  momentum OutR = Var.LoopMom[OUTR];
  momentum K0 = Var.LoopMom[4];
  double Factor = 1.0 / pow(2.0 * PI, D);
  double dTau, GWeight, GWeightInBox;
  verWeight Weight, CWeight, RefWeight, TestWeight;

  // one loop T and TC diagram
  momentum K1 = OutL + K0 - InL;
  auto LVer = Prop.Interaction(InL, OutL, K1, K0, false);
  auto RVer = Prop.Interaction(K0, K1, InR, OutR, false);
  dTau = Var.Tau[1] - Var.Tau[0];
  GWeight = Prop.Green(dTau, K0, UP, 0) * Prop.Green(-dTau, K1, UP, 0);
  GWeightInBox = Prop.Green(dTau, K0, UP, 0) * Prop.Green(-dTau, K0, UP, 0);

  Weight[DIR] = LVer[DIR] * RVer[DIR] * SPIN + LVer[EX] * RVer[DIR] +
                LVer[DIR] * RVer[EX];
  Weight[EX] = LVer[EX] * RVer[EX];
  Weight *= GWeight;
  Weight *= Factor * SymFactor[T];

  CWeight[DIR] = GWeightInBox * LVer[DIR] * RVer[DIR] * Para.Lambda /
                 (8.0 * PI) / Para.Nf * SPIN;
  CWeight[EX] = 0.0;
  // cout << CWeight[DIR] << endl;
  CWeight *= Factor * SymFactor[T];

  TestWeight = Weight - CWeight;

  Weight = _GetWeight(1, {T, TC});

  if (abs(Weight[DIR] - TestWeight[DIR]) > 1.0e-10 ||
      abs(Weight[EX] - TestWeight[EX]) > 1.0e-10) {
    cout << fmt::format("G0: {}", Prop.Green(dTau, K0, UP, 0)) << endl;
    cout << fmt::format("GT: {}", Prop.Green(-dTau, K1, UP, 0)) << endl;
    cout << fmt::format("G: {}", GWeight) << endl;
    cout << fmt::format("Gbox: {}", GWeightInBox) << endl;
    // cout << fmt::format(
    //             "GInBox: {}, sep: {}, {}",
    //             GWeightInBox * Factor * Para.Lambda / (8.0 * PI) / Para.Nf,
    //             Prop.Green(dTau, K0, UP, 0), Prop.Green(-dTau, K0, UP, 0))
    //      << endl;
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
  TestWeight *= Factor * SymFactor[S] * cos(2.0 * PI / Para.Beta * dTau * 2.0);
  // TestWeight *= Factor * SymFactor[S];

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

verWeight vertex4::_GetWeight(int LoopNum, std::vector<channel> Channel) {
  verWeight Weight;
  // array<momentum *, 4> ExtLegK = {&Var.LoopMom[0], &Var.LoopMom[1],
  //                                 &Var.LoopMom[2], &Var.LoopMom[3]};
  vertex4 Ver4;
  Ver4.Build(0, // level
             1, // loopNum
             4, // loop index of the first internal K
             0, // tau index of the InTL leg
             Channel, RIGHT);

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
