#include "global.h"
#include "utility/abort.h"
#include "utility/fmt/format.h"
#include "utility/timer.h"
#include "utility/vector.h"
#include "weight.h"
#include <array>
#include <iostream>
#include <stack>
#include <string>

using namespace diag;
using namespace std;
using namespace dse;

void weight::Test(int LoopNum, diagram Diagram) {

  LOG_INFO("Testing T channel one loop ... ");
  vector<channel> Chan = {T, TC};
  auto Weight = _GetWeight(1, Chan);
  momentum InL = Var.LoopMom[INL];
  momentum OutL = Var.LoopMom[OUTL];
  momentum InR = Var.LoopMom[INR];
  momentum OutR = Var.LoopMom[OUTR];
  momentum K0 = Var.LoopMom[4];

  // one loop T and TC diagram
  momentum K1 = OutL + K0 - InL;
  array<momentum *, 4> LLegK = {&InL, &OutL, &K0, &K1};
  auto LVer = VerQTheta.Interaction(LLegK, 0, false, false);
  array<momentum *, 4> RLegK = {&K0, &K1, &InR, &OutR};
  auto RVer = VerQTheta.Interaction(RLegK, 0, false, false);
  double dTau = Var.Tau[1] - Var.Tau[0];
  double GWeight =
      Fermi.Green(dTau, K0, UP, 0, 0.0) * Fermi.Green(-dTau, K1, UP, 0, 0.0);
  double GWeightInBox =
      Fermi.Green(dTau, K0, UP, 0, 0.0) * Fermi.Green(-dTau, K0, UP, 0, 0.0);

  double Factor = 1.0 / pow(2.0 * PI, D);
  ver::weightMatrix TestWeight;
  TestWeight.SetZero();
  TestWeight[DIR] = LVer[DIR] * RVer[DIR] * SPIN + LVer[EX] * RVer[DIR] +
                    LVer[DIR] * RVer[EX];
  TestWeight[EX] = LVer[EX] * RVer[EX];
  TestWeight *= GWeight;
  momentum TranQ = OutL - InL;
  TestWeight[DIR] -=
      GWeightInBox *
      pow(8.0 * PI / (TranQ.squaredNorm() + Para.Mass2 + Para.Lambda), 2) *
      Para.Lambda / (8.0 * PI) / Para.Nf * GWeightInBox * SPIN;
  TestWeight *= Factor * SymFactor[T];
  if (abs(Weight[DIR] - TestWeight[DIR]) > 1.0e-10 ||
      abs(Weight[EX] - TestWeight[EX]) > 1.0e-10) {
    // cout << fmt::format("G: {}", GWeight * Factor) << endl;
    // cout << fmt::format("LVer: {}, {}", LVer[DIR], LVer[EX]) << endl;
    // cout << fmt::format("RVer: {}, {}", RVer[DIR], RVer[EX]) << endl;
    ABORT(fmt::format("Weight Calculation has a bug! ({},{}) vs ({}, {})",
                      Weight[DIR], Weight[EX], TestWeight[DIR],
                      TestWeight[EX]));
  }
  return;
}

ver::weightMatrix weight::_GetWeight(int LoopNum, vector<channel> Channel) {
  ver::weightMatrix Weight;
  array<momentum *, 4> ExtLegK = {&Var.LoopMom[0], &Var.LoopMom[1],
                                  &Var.LoopMom[2], &Var.LoopMom[3]};
  ver4 Ver4 = VerDiag.Vertex(0, // level
                             1, // loopNum
                             4, // loop index of the first internal K
                             0, // tau index of the InTL leg
                             Channel, RIGHT, false);
  VerDiag.ResetMomMap(Ver4, ExtLegK);
  Vertex4(Ver4, true);
  for (auto &w : ChanWeight)
    // collapse all channel to I
    ChanWeight[0] += w;
  // cout << ChanWeight[0].Abs() << endl;
  Weight = ChanWeight[0];
  return Weight;
}
